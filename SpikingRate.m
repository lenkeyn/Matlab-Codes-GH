function sData = SpikingRate(sData)

signal = sData.imdata.roiSignals(2).deconv;
fps = sData.imdata.meta.fps;
nROIs = size(signal,1);
nSamples = size(signal,2);

%%% ESTIMATE 'SPIKING RATE', or cell activity (arbitrary 'firing rate') using the deconved signal
% By integrating the deconvolved signal (it looks weird spiky reporting each activation) within a time window and setting a treshold what to consider activation, we can calculate some kind of arbitrary spike rate. At the end we apply gaussian smooting.(?) Credit to Eivind Hennestad for most of this, and the remaining credit is for Andreas...

% Set parmateres:
CellFactor = 0.25; %multiplying factor for different cell tyes, PC = 0.5 or 1; VIP = 2;
cumSamples = ceil(CellFactor*sData.imdata.meta.fps); % how many samples to use for 'memory' (cumulate signal) (31 frames = 1 sec). I set 0.5 sec (15 samples) for pyramidal cells, because rise of Ca-transient is usually within 0.5 sec. VIP: 2 sec
TestedSampleIndex = round(cumSamples/2);  % the tested sample should be in the middle of the activity array
%spikethresholdArray = max(sData.imdata.roiSignals(2).deconv,[],2); % treshold for consider an event as 'firing' (to discard small and fast fluorescent changes which does not seem real activation although the deconvolution method found them). Eivind uses 0.24, Andreas 0.11. I set 0.05, I found it reasonably good (maybe not detecting all single action potentials, but most of it does). Increase it to 0.24 might cause loosing 'few action potential' events, decreasing it to 0.01 might cause detection of a few noise deflections in the signal detected by deconvolution.
spikethreshold = 0.05; %0.25 , 0.1, 0.05
estimated_firing_rate = zeros(nROIs,nSamples); % create empty array for data

% calculating spike rate
for i = 1:1:nROIs
    %spikethreshold = spikethresholdArray(i)*0.2; 
    SignalROIdeconv = sData.imdata.roiSignals(2).deconv(i,1:nSamples); % temporary array for calculating firing rate
    if isnan(SignalROIdeconv(i))
        continue
    end
    cumSum = 0; % The sum of integrating the spike rate
    SmoothingWindow = zeros(1,cumSamples); % Temporary array , activity around sample. I want to reset the sum if nothing happened for the last second. Could maybe be shorter period.
    SignalROIspikes = zeros(size(SignalROIdeconv)); % Create vector to put spikes in
    SignalROIspikes(1:TestedSampleIndex-1) = SignalROIdeconv(1:TestedSampleIndex-1); % keep the deconvolved signal in the beginning and at the end
    SignalROIspikes(end-TestedSampleIndex:end) = SignalROIdeconv(1:TestedSampleIndex-1);
    
    for j = 1:1:30 %length(SignalROIdeconv) % Run through the deconvolved signal
        SmoothingWindow = horzcat(SmoothingWindow(2:end),SignalROIdeconv(j)); % Put current sample at the end of "memory" vector
        if sum(SmoothingWindow) == 0  % if there are only zeros, reset the cumSum.
            cumSum = 0;
            continue
        elseif SignalROIdeconv(j) > 0 % if there is a value in the deconvolved signal, keep it
            SignalROIspikes(j) = SignalROIdeconv(j);
            continue
        else  % if there are zeros between two or more values, use the mean of the values
            NonZeroIndices = SmoothingWindow~=0;
            NonZeroValues = SmoothingWindow(NonZeroIndices);
            cumSum = mean(NonZeroValues);
            if  cumSum/spikethreshold >= 1 % Check if sum is over threshold. Count spikes and reset.
                SignalROIspikes(j) = cumSum; %  SignalROIspikes(j) = floor(cumSum/spikethreshold); 
            else
                cumSum = 0;
            end 
        end
        
        %{
        if floor(cumSum/spikethreshold) >= 1 % Check if sum is over threshold. Count spikes and reset.
            SignalTemp4(j) = floor(cumSum/spikethreshold); %  SignalTemp4(j) = floor(cumSum/spikethreshold); 
            cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold; % cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold;
        end
        %}
    end
    estimated_firing_rate(i,:) = SignalROIspikes;
    %estimated_firing_rate(i,:) = smoothdata(SignalROI2,'gaussian',round(fps/2));  %PC: round(FrameRate/2) % old:smoothdata(SignalTemp4./(1/FrameRate),'gaussian',round(FrameRate));smoothdata(SignalTemp4./(1/FrameRate),'gaussian',round(FrameRate));
end
sData.imdata.roiSignals(2).SpikingRate = single(estimated_firing_rate);

end
