function behav = OptoTrialSorting2(behav,sData,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm)
% TDMSDATA has been changed to sData.daqdata 20190221
% Sensitivity = 100; %% If there was stimulation in the (1/Sensitivity)percentage during the trial it should be considered as opto-on trial. 1/100 high sens, 1/10 low
savePath = strcat(sData.sessionInfo.savePath,'\Behavior');

% Set time limit for optical stimulation. If optical stimulation lasts longer,discard that trial
%OptoStimLimit = 9000; % millisec
%DiscardCmBeginningCm = 10; %do not count the first e.g.20 cms (in bins), because many time stimulation from previous trial ends
% set a treshold to detect optical stimulation: 50% of peak stimulation
PhotoStimLimit = (max(sData.daqdata.optoSignal) - min(sData.daqdata.optoSignal))/100; % size of photostim/100
% generate a binned opto signal, zero: no stimulation, one: stimulation
OptoOnSignal = zeros(size(sData.daqdata.optoSignal,1),1);
OptoOnSignal(sData.daqdata.optoSignal > PhotoStimLimit) = 1; % if there is opto stimulation set value to 1, otherwise zero
% generate downsampled opto signal
behav.Opto.OptoOnSignalDS = OptoOnSignal(behav.FrameStartIndex);

TRNu = size(behav.LickBinMatrixExtended,1);
% collect the trial indices into four arrays: opto-on, opto-off, after-opto, failed opto trials
behav.Opto.OptoOnTrials = NaN(TRNu,1); 
behav.Opto.OptoOffTrials = NaN(TRNu,1); 
behav.Opto.AfterOptoTrials = NaN(TRNu,1);
behav.Opto.FailedOptoTrials = NaN(TRNu,1);
StartingBin = ceil(DiscardCmBeginningCm/behav.BinSize)+1; % do not count the first e.g.20 cms (in bins), because many time stimulation from previous trial ends
counter = 0;
for i = 1:1:TRNu
    TrialOpto = behav.Opto.OptoOnSignalDS(behav.EnterIntoBinSampleInd(i,StartingBin):behav.EnterIntoBinSampleInd(i+1,1)-1); 
    if sum(TrialOpto)< numel(TrialOpto)/OptoSensitivity  % If there was Opto in the (1/Sensitivity) of the trial it should be considered as Opto-on trial
        behav.Opto.OptoOffTrials(i) = 1;
    else
        % discard trial if opto stimulation was longer then 9 seconds (9000 ms)
        TrialOptoFullTrial = behav.Opto.OptoOnSignalDS(behav.EnterIntoBinSampleInd(i,1):behav.EnterIntoBinSampleInd(i+1,1)-1); 
        TrialOptoFullTrialStartIndex = find(TrialOptoFullTrial~=0, 1 );
        TrialOptoFullTrialEndIndex = find(TrialOptoFullTrial~=0, 1, 'last');
        OptoStimDuration = (TrialOptoFullTrialEndIndex - TrialOptoFullTrialStartIndex) * (1000/sData.behavior.meta.imagingSamplingRate); % optical stimulation in ms
        if OptoStimDuration > OptoStimLimitMs
            behav.Opto.OptoOnTrials(i) = NaN; % too long optical stimulation
            behav.Opto.FailedOptoTrials(i) = 1;
            counter = counter +1;
            continue
        end
        behav.Opto.OptoOnTrials(i) = 1; % correct opto-on trial
    end
end
msgbox(sprintf('Opto stimulation was too long in trials: %g', counter));
       
% search the indices of after-opto trials
OptoOnOffChange = behav.Opto.OptoOnTrials;
OptoOnOffChange(isnan(OptoOnOffChange)) = 0;
OptoOnOff = find(diff(OptoOnOffChange)== -1);
AfterOptoTrials = OptoOnOff+1;
behav.Opto.AfterOptoTrials(AfterOptoTrials) = 1;

behav.Opto.OptoOffTrials(behav.Opto.AfterOptoTrials == 1) = 0;

behav.Opto.OptoOffTrialsIndices = find(behav.Opto.OptoOffTrials == 1);
behav.Opto.OptoOnTrialsIndices = find(behav.Opto.OptoOnTrials == 1);
behav.Opto.AfterOptoTrialsIndices = find(behav.Opto.AfterOptoTrials == 1);
behav.Opto.FailedOptoTrialsIndices = find(behav.Opto.FailedOptoTrials == 1);

%%% Sort the licking data
behav.Opto.LickMatrix_OptoOff = behav.LickCmMatrixExtended(behav.Opto.OptoOffTrialsIndices,:);
behav.Opto.LickMatrix_OptoOn = behav.LickCmMatrixExtended(behav.Opto.OptoOnTrialsIndices,:);
behav.Opto.LickMatrix_AfterOpto = behav.LickCmMatrixExtended(behav.Opto.AfterOptoTrialsIndices,:);

behav.Opto.MeanLick_OptoOff = nanmean(behav.Opto.LickMatrix_OptoOff,1);
behav.Opto.MeanLick_OptoOn = nanmean(behav.Opto.LickMatrix_OptoOn,1);
behav.Opto.MeanLick_OptoAfter = nanmean(behav.Opto.LickMatrix_AfterOpto,1);

%%%% Sort speed data
behav.Opto.VeloMatrix_OptoOff = behav.VeloBinMatrixExtended(behav.Opto.OptoOffTrialsIndices,:);
behav.Opto.VeloMatrix_OptoOn = behav.VeloBinMatrixExtended(behav.Opto.OptoOnTrialsIndices,:);
behav.Opto.VeloMatrix_AfterOpto = behav.VeloBinMatrixExtended(behav.Opto.AfterOptoTrialsIndices,:);

behav.Opto.MeanVelo_OptoOff = mean(behav.Opto.VeloMatrix_OptoOff,1);
behav.Opto.MeanVelo_OptoOn = mean(behav.Opto.VeloMatrix_OptoOn,1);
behav.Opto.MeanVelo_OptoAfter = mean(behav.Opto.VeloMatrix_AfterOpto,1);

%%% HitRate
PosNaN = find(isnan(behav.Opto.HitTrialsArray));
behav.Opto.ActiveHitRate_OptoOff = sum(behav.Opto.HitTrialsArray(behav.Opto.OptoOffTrialsIndices(behav.Opto.OptoOffTrialsIndices<PosNaN),1)/nnz(~isnan(behav.Opto.HitTrialsArray(behav.Opto.OptoOffTrialsIndices,1))));
behav.Opto.ActiveHitRate_OptoOn = sum(behav.Opto.HitTrialsArray(behav.Opto.OptoOnTrialsIndices(behav.Opto.OptoOnTrialsIndices<PosNaN),1)/nnz(~isnan(behav.Opto.HitTrialsArray(behav.Opto.OptoOnTrialsIndices,1))));
behav.Opto.ActiveHitRate_OptoAfter = sum(behav.Opto.HitTrialsArray(behav.Opto.AfterOptoTrialsIndices(behav.Opto.AfterOptoTrialsIndices<PosNaN),1)/nnz(~isnan(behav.Opto.HitTrialsArray(behav.Opto.AfterOptoTrialsIndices,1))));

%%% Binning Opto-stimulus signal to check in which position was opto stimulus
behav.Opto.OptoStimOnMatrix = NaN(TRNu,behav.BinNu);
% detect the optical stim during each bin 
for i = 1:1:TRNu  % rows are trials
    for j = 1:1:behav.BinNu  % columns (distance bins)  
        BinnedOpto = mean(behav.Opto.OptoOnSignalDS(behav.EnterIntoBinSampleInd(i,j):behav.LeaveBinSampleInd(i,j)));
        if BinnedOpto > 0.001 % I set the treshold low to detect all stimulations, sometimes there are brief false positive spikes in the signal, therefore I need to treshold it
            behav.Opto.OptoStimOnMatrix(i,j) = 1;  
        else
            behav.Opto.OptoStimOnMatrix(i,j) = 0;
        end
    end
end


%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:TRNu,behav.Opto.OptoStimOnMatrix) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
mymap = [0.8 0.8 0.8
         0.9 0.2 0.2]; % set colormap colors, gray and red
colormap(mymap);
xlabel('Position on wheel (cm)');
ax = gca; ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(strcat(behav.FileID,'-optical-stim'));
FileName = strcat('OptoStimHeatBin-',behav.FileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

StimMax = max(mean(behav.Opto.OptoStimOnMatrix));
behav.Opto.optoStimStart = min(find(mean(behav.Opto.OptoStimOnMatrix)>= StimMax/2,1)*behav.BinSize);
behav.Opto.optoStimEnd = max(find(mean(behav.Opto.OptoStimOnMatrix)>= StimMax/2,1,'last')*behav.BinSize);

end