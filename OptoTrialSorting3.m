function Opto = OptoTrialSorting3(sData,OptoSensitivity,OptoStimLimitMs,DiscardBeginningCm)
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
Opto.OptoOnSignalDS = OptoOnSignal(sData.behavior.details.frameStartIndices);

TRNu = sData.behavior.WheelLapImagingExtended; %LickBinMatrixExtended
% collect the trial indices into four arrays: opto-on, opto-off, after-opto, failed opto trials
Opto.OptoOnTrials = NaN(TRNu,1); 
Opto.OptoOffTrials = NaN(TRNu,1); 
Opto.AfterOptoTrials = NaN(TRNu,1);
Opto.FailedOptoTrials = NaN(TRNu,1);
StartingBin = ceil(DiscardBeginningCm/sData.behavior.meta.binSize)+1; % do not count the first e.g.20 cms (in bins), because many time stimulation from previous trial ends
counter = 0;
for i = 1:1:TRNu
    TrialOpto = Opto.OptoOnSignalDS(sData.behavior.binning.enterIntoBinIndex(i,StartingBin):sData.behavior.binning.enterIntoBinIndex(i+1,1)-1); 
    if sum(TrialOpto)< numel(TrialOpto)/OptoSensitivity  % If there was Opto in the (1/Sensitivity) of the trial it should be considered as Opto-on trial
        %Opto.OptoOffTrials(i) = 1;
    else
        % discard trial if opto stimulation was longer then 9/15 seconds (9000 ms)
        TrialOptoFullTrial = Opto.OptoOnSignalDS(sData.behavior.binning.enterIntoBinIndex(i,1):sData.behavior.binning.enterIntoBinIndex(i+1,1)-1); 
        TrialOptoFullTrialStartIndex = find(TrialOptoFullTrial~=0, 1 );
        TrialOptoFullTrialEndIndex = find(TrialOptoFullTrial~=0, 1, 'last');
        OptoStimDuration = (TrialOptoFullTrialEndIndex - TrialOptoFullTrialStartIndex) * (1000/sData.behavior.meta.imagingSamplingRate); % optical stimulation in ms
        if OptoStimDuration > OptoStimLimitMs
            Opto.OptoOnTrials(i) = NaN; % too long optical stimulation
            Opto.FailedOptoTrials(i) = 1;
            counter = counter +1;
            continue
        end
        Opto.OptoOnTrials(i) = 1; % correct opto-on trial
    end
end
msgbox(sprintf('Opto stimulation was too long in trials: %g', counter));
       
% search the indices of after-opto trials
OptoOnOffChange = Opto.OptoOnTrials;
OptoOnOffChange(isnan(OptoOnOffChange)) = 0;
OptoOnOff = find(diff(OptoOnOffChange)== -1);
AfterOptoTrials = OptoOnOff+1;
Opto.AfterOptoTrials(AfterOptoTrials) = 1;
Opto.AfterOptoTrials(Opto.FailedOptoTrials==1) = NaN;
Opto.AfterOptoTrials(find(Opto.FailedOptoTrials==1)+1) = NaN;

Opto.OptoOffTrials(:,1) = 1;
Opto.OptoOffTrials(Opto.AfterOptoTrials == 1) = 0;
Opto.OptoOffTrials(Opto.FailedOptoTrials==1) = NaN;
Opto.OptoOffTrials(Opto.OptoOnTrials == 1) = NaN;

Opto.OptoOffTrialsIndices = find(Opto.OptoOffTrials == 1);
Opto.OptoOnTrialsIndices = find(Opto.OptoOnTrials == 1);
Opto.AfterOptoTrialsIndices = find(Opto.AfterOptoTrials == 1);
Opto.FailedOptoTrialsIndices = find(Opto.FailedOptoTrials == 1);

%%% Sort the licking data
Opto.LickMatrix_OptoOff = sData.behavior.binning.lickPerCmBinned(Opto.OptoOffTrialsIndices,:);
Opto.LickMatrix_OptoOn = sData.behavior.binning.lickPerCmBinned(Opto.OptoOnTrialsIndices,:);
Opto.LickMatrix_AfterOpto = sData.behavior.binning.lickPerCmBinned(Opto.AfterOptoTrialsIndices,:);

Opto.MeanLick_OptoOff = nanmean(Opto.LickMatrix_OptoOff,1);
Opto.MeanLick_OptoOn = nanmean(Opto.LickMatrix_OptoOn,1);
Opto.MeanLick_OptoAfter = nanmean(Opto.LickMatrix_AfterOpto,1);

%%%% Sort speed data
Opto.VeloMatrix_OptoOff = sData.behavior.binning.veloBinned(Opto.OptoOffTrialsIndices,:);
Opto.VeloMatrix_OptoOn = sData.behavior.binning.veloBinned(Opto.OptoOnTrialsIndices,:);
Opto.VeloMatrix_AfterOpto = sData.behavior.binning.veloBinned(Opto.AfterOptoTrialsIndices,:);

Opto.MeanVelo_OptoOff = mean(Opto.VeloMatrix_OptoOff,1);
Opto.MeanVelo_OptoOn = mean(Opto.VeloMatrix_OptoOn,1);
Opto.MeanVelo_OptoAfter = mean(Opto.VeloMatrix_AfterOpto,1);

%%% HitRate
PosNaN = find(isnan(sData.behavior.details.hitTrials));
Opto.ActiveHitRate_OptoOff = sum(sData.behavior.details.hitTrials(Opto.OptoOffTrialsIndices(Opto.OptoOffTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(Opto.OptoOffTrialsIndices,1))));
Opto.ActiveHitRate_OptoOn = sum(sData.behavior.details.hitTrials(Opto.OptoOnTrialsIndices(Opto.OptoOnTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(Opto.OptoOnTrialsIndices,1))));
Opto.ActiveHitRate_OptoAfter = sum(sData.behavior.details.hitTrials(Opto.AfterOptoTrialsIndices(Opto.AfterOptoTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(Opto.AfterOptoTrialsIndices,1))));

%%% Binning Opto-stimulus signal to check in which position was opto stimulus
Opto.OptoStimOnMatrix = NaN(TRNu,sData.behavior.meta.nBins);
% detect the optical stim during each bin 
for i = 1:1:TRNu  % rows are trials
    for j = 1:1:sData.behavior.meta.nBins  % columns (distance bins)  
        BinnedOpto = mean(Opto.OptoOnSignalDS(sData.behavior.binning.enterIntoBinIndex(i,j):sData.behavior.binning.leaveBinIndex(i,j)));
        if BinnedOpto > 0.001 % I set the treshold low to detect all stimulations, sometimes there are brief false positive spikes in the signal, therefore I need to treshold it
            Opto.OptoStimOnMatrix(i,j) = 1;  
        else
            Opto.OptoStimOnMatrix(i,j) = 0;
        end
    end
end


%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:TRNu,Opto.OptoStimOnMatrix) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
mymap = [0.8 0.8 0.8
         0.9 0.2 0.2]; % set colormap colors, gray and red
colormap(mymap);
xlabel('Position on wheel (cm)');
ax = gca; ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(strcat(sData.sessionInfo.fileID,'-optical-stim'));
FileName = strcat('OptoStimHeatBin-',sData.sessionInfo.fileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

StimMax = max(mean(Opto.OptoStimOnMatrix));
Opto.optoStimStart = min(find(mean(Opto.OptoStimOnMatrix)>= StimMax/2,1)*sData.behavior.meta.binSize);
Opto.optoStimEnd = max(find(mean(Opto.OptoStimOnMatrix)>= StimMax/2,1,'last')*sData.behavior.meta.binSize);

end