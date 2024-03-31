function sData = cumSpikes(sData)

signal = sData.imdata.roiSignals(2).deconv;
fps = sData.imdata.meta.fps;
nROIs = size(signal,1);
nSamples = size(signal,2);

%%% ESTIMATE 'SPIKING RATE', or cell activity (arbitrary 'firing rate') using the deconved signal
% By integrating the deconvolved signal (it looks weird spiky reporting each activation) within a time window and setting a treshold what to consider activation, we can calculate some kind of arbitrary spike rate. At the end we apply gaussian smooting.(?) Credit to Eivind Hennestad for most of this, and the remaining credit is for Andreas...

% Set parmateres:
cumTime = 500; % time in ms to cumulate signal (for PC: 200-300 ms, for VIP more)
cumSamples = floor(cumTime/fps); % how many samples to use for 'memory' (cumulate signal) (31 frames = 1 sec). I set 0.5 sec (15 samples) for pyramidal cells, because rise of Ca-transient is usually within 0.5 sec. VIP: 2 sec
TestedSampleIndex = round(cumSamples/2);  % the tested sample should be in the middle of the activity array
%spikethresholdArray = max(sData.imdata.roiSignals(2).deconv,[],2); % treshold for consider an event as 'firing' (to discard small and fast fluorescent changes which does not seem real activation although the deconvolution method found them). Eivind uses 0.24, Andreas 0.11. I set 0.05, I found it reasonably good (maybe not detecting all single action potentials, but most of it does). Increase it to 0.24 might cause loosing 'few action potential' events, decreasing it to 0.01 might cause detection of a few noise deflections in the signal detected by deconvolution.
spikethreshold = 0.05; %0.25 , 0.1, 0.05
estimated_firing_rate = zeros(nROIs,nSamples); % create empty array for data

% calculating spike rate
for i = 1:1:nROIs
    %spikethreshold = spikethresholdArray(i)*0.2; % for aoutomatic estimation
    SignalROIdeconv = sData.imdata.roiSignals(2).deconv(i,1:nSamples); % temporary array for calculating firing rate
    if isnan(SignalROIdeconv(i))
        continue
    end
    
    % first round:
    SignalROIspikes1 = zeros(size(SignalROIdeconv)); % Create vector to put spikes in for a ROI
    SignalROIspikes1(1:TestedSampleIndex-1) = SignalROIdeconv(1:TestedSampleIndex-1); % keep the deconvolved signal in the beginning and at the end
    SignalROIspikes1(end-TestedSampleIndex:end) = SignalROIdeconv(end-TestedSampleIndex:end);
    for j = TestedSampleIndex:1:length(SignalROIdeconv)-TestedSampleIndex-1 % Run through the deconvolved signal
        SmoothingWindow = SignalROIdeconv(j-TestedSampleIndex+1:j+TestedSampleIndex-1); % Put current sample in the middle 
        if SignalROIdeconv(j) > 0 % if there is a value in the deconvolved signal, keep it
            SignalROIspikes1(j) = SignalROIdeconv(j);
            continue
        elseif sum(SmoothingWindow(1:TestedSampleIndex)) == 0 || sum(SmoothingWindow(TestedSampleIndex:end)) == 0  % if there are only zeros before or after the sample, keep the value at zero.
            continue
        else  % if the value is zero between two or more values, use the mean of the values
            NonZeroIndices = SmoothingWindow~=0;
            NonZeroValues = SmoothingWindow(NonZeroIndices);
            cumSum = mean(NonZeroValues);
            if  cumSum/spikethreshold >= 1 % Check if sum is over threshold. If yes update value, if not keep zero.
                SignalROIspikes1(j) = cumSum;  
            end 
        end
    end
    
    %second round:
    SignalROIspikes2 = zeros(size(SignalROIdeconv)); % Create vector to put spikes in for a ROI
    SignalROIspikes2(1:TestedSampleIndex-1) = SignalROIdeconv(1:TestedSampleIndex-1); % keep the deconvolved signal in the beginning and at the end
    SignalROIspikes2(end-TestedSampleIndex:end) = SignalROIdeconv(end-TestedSampleIndex:end);
    for j = TestedSampleIndex:1:length(SignalROIdeconv)-TestedSampleIndex-1 % Run through the deconvolved signal
        SmoothingWindow = SignalROIspikes1(j-TestedSampleIndex+1:j+TestedSampleIndex-1); % Put current sample in the middle 
        if SignalROIspikes1(j) > 0 % if there is a value in the deconvolved signal, keep it
            SignalROIspikes2(j) = SignalROIspikes1(j);
            continue
        elseif sum(SmoothingWindow(1:TestedSampleIndex)) == 0 || sum(SmoothingWindow(TestedSampleIndex:end)) == 0  % if there are only zeros before or after the sample, keep the value at zero.
            continue
        else  % if the value is zero between two or more values, use the mean of the values
            NonZeroIndices = SmoothingWindow~=0;
            NonZeroValues = SmoothingWindow(NonZeroIndices);
            cumSum = mean(NonZeroValues);
            if  cumSum/spikethreshold >= 1 % Check if sum is over threshold. If yes update value, if not keep zero.
                SignalROIspikes2(j) = cumSum;  
            end 
        end
    end
    % collecting firing rate data
    estimated_firing_rate(i,:) = SignalROIspikes2;
    
end

sData.imdata.roiSignals(2).cumSpikes = single(estimated_firing_rate);

end
