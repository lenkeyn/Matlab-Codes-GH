function sData = SlowRemoved(sData,Window1,Percentile1)

%%% Remove slow time-scale changes in fluoresence and baseline subtraction (5 percentile), baseline subtraction is not perfect (always above zero)
%Inspired Andreas by Dombeck et al (2010, Nat Neuro), code written by Andreas. Calculate basline in every 15 sec (in a moveing window) and subtract from data
RawNpilSubt = sData.imdata.roiSignals(2).roif - sData.imdata.roiSignals(2).npilf; 

FrameRate = 31;
nROIs = sData.imdata.nROIs;
nSamples = sData.imdata.nSamples;
%Window1 = 10; % seconds, can be changed... , 10-15s, start: 10s
%Percentile1 = 20; % I used 20%, but the transients became negative in baseline
Sampl = ceil(Window1*FrameRate); % samples to be used
SignalTemp1 = zeros(nROIs,nSamples); % temporary array for signal calculation
MeanBaseline = zeros(nROIs,1); % mean baseline for each ROI recording for calculate dFF
CollectBaseline = zeros(1,nSamples); % temporary collection of baseline value for a ROI
for i = 1:1:nROIs 
    % duplicate the beginning and end of the signal, and concatenate the first (and last) 15s of the signal to the beginning (and end) of the original signal
    SignalTemp2 = [RawNpilSubt(i,1:Sampl),RawNpilSubt(i,:),RawNpilSubt(i,(end-Sampl+1):end)]; % concatenate the duplicates and the original signal
    % Use a window of Window seconds around each data point to obtain X th percentile and subtract this from original signal to make the basline flat
    for j = Sampl:1:(nSamples+Sampl-1)
        Signal_window = SignalTemp2(1,(j-round(Sampl/2)):(j+round(Sampl/2))); % collect the data witihn the actual window to calculate baseline for datapoint 
        SignalTemp1(i,(j-Sampl+1)) = SignalTemp2(j) - prctile(Signal_window,Percentile1); % calculate X percentile of data and subtract
        CollectBaseline(1,(j-Sampl+1)) = prctile(Signal_window,Percentile1);
    end
    MeanBaseline(i,1) = mean(CollectBaseline(1,:)); % mean baseline for the whole ROI session, used for dF/F division
end
ROIsignals_raw_slow_removed = SignalTemp1; %

% Baseline already subtracted, generate dFF (devide with baseline)
dFF_slowRemoved = ROIsignals_raw_slow_removed ./ abs(MeanBaseline); % calculate dFF (baseline was subtracted in previous session)
sData.imdata.roiSignals(2).dff_slowRemoved = single(dFF_slowRemoved);

end


%{
roi = 11;
figure (1)
plot(movmean(sData.imdata.roiSignals(2).dff(roi,1:5000),20)); hold on
plot(movmean(dFF_slowRemoved3p20(roi,1:5000)+2,20)); hold on
plot(movmean(dFF_slowRemoved10p20(roi,1:5000)+4,20)); hold on
plot(movmean(dFF_slowRemoved30p20(roi,1:5000)+6,20)); hold on
%}


% save temp
%save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
