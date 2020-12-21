function LVDATA = OptoTrialSorting2(LVDATA,sData,Sensitivity)
% TDMSDATA has been changed to sData.daqdata 20190221
%Sensitivity = 100; %% If there was light in the (1/Sensitivity)percentage during the trial it should be considered as light-on trial. 1/100 high sens, 1/10 low
savePath = strcat(sData.sessionInfo.savePath,'\Behavior');

PhotoStimLimit = (max(sData.daqdata.optoSignal) - min(sData.daqdata.optoSignal))/2; % size of photostim
LightOnSignal = zeros(size(sData.daqdata.optoSignal,1),1);
LightOnSignal(sData.daqdata.optoSignal > PhotoStimLimit) = 1;
LVDATA.Opto.LightOnSignalDS = LightOnSignal(LVDATA.FrameStartIndex);

TRNu = size(LVDATA.LickBinMatrix,1);
LVDATA.Opto.LightOnTrials = NaN(TRNu,1); 
LVDATA.Opto.LightOffTrials = NaN(TRNu,1); 
LVDATA.Opto.AfterLightTrials = NaN(TRNu,1);
for i = 1:1:TRNu
    %TrialLight = NaN(max(diff(LVDATA.EnterIntoBinSampleInd(:,1)))+1,1); % temporary array for a trial 
    TrialLight = LVDATA.Opto.LightOnSignalDS(LVDATA.EnterIntoBinSampleInd(i,4):LVDATA.EnterIntoBinSampleInd(i+1,1)-1); % do not count the first 3 bins, because many time stimulation from previous trial ends in the next
    if sum(TrialLight)< numel(TrialLight)/Sensitivity  % If there was light in the (1/Sensitivity) of the trial it should be considered as light-on trial
        LVDATA.Opto.LightOffTrials(i) = 1;
    else
        if numel(TrialLight)> 450 % if trial lasats longer then 15 sec, discard it
            continue 
        end
        LVDATA.Opto.LightOnTrials(i) = 1;
    end
end
LightOnOffChange = LVDATA.Opto.LightOnTrials;
LightOnOffChange(isnan(LightOnOffChange)) = 0;
LightOnOff = find(diff(LightOnOffChange)== -1);
AfterLightTrials = LightOnOff+1;
LVDATA.Opto.AfterLightTrials(AfterLightTrials) = 1;
LVDATA.Opto.LightOffTrials(LVDATA.Opto.AfterLightTrials == 1) = 0;

LightOffTrials = find(LVDATA.Opto.LightOffTrials == 1);
LightOnTrials = find(LVDATA.Opto.LightOnTrials == 1);
AfterLightTrials = find(LVDATA.Opto.AfterLightTrials == 1);

%%% NORI EXTRA

%%%Lick
LVDATA.Opto.LickMatrix_LightOff = LVDATA.LickCmMatrix(LightOffTrials,:);
LVDATA.Opto.LickMatrix_LightOn = LVDATA.LickCmMatrix(LightOnTrials,:);
LVDATA.Opto.LickMatrix_AfterLight = LVDATA.LickCmMatrix(AfterLightTrials,:);

LVDATA.Opto.MeanLick_LightOff = nanmean(LVDATA.Opto.LickMatrix_LightOff);
LVDATA.Opto.MeanLick_LightOn = nanmean(LVDATA.Opto.LickMatrix_LightOn);
LVDATA.Opto.MeanLick_LightAfter = nanmean(LVDATA.Opto.LickMatrix_AfterLight);


%%%% SPEED
LVDATA.Opto.VeloMatrix_LightOff = LVDATA.VeloBinMatrix(LightOffTrials,:);
LVDATA.Opto.VeloMatrix_LightOn = LVDATA.VeloBinMatrix(LightOnTrials,:);
LVDATA.Opto.VeloMatrix_AfterLight = LVDATA.VeloBinMatrix(AfterLightTrials,:);

LVDATA.Opto.MeanVelo_LightOff = mean(LVDATA.Opto.VeloMatrix_LightOff);
LVDATA.Opto.MeanVelo_LightOn = mean(LVDATA.Opto.VeloMatrix_LightOn);
LVDATA.Opto.MeanVelo_LightAfter = mean(LVDATA.Opto.VeloMatrix_AfterLight);

%%% HitRate
PosNaN = find(isnan(LVDATA.Opto.HitTrialsArray));
LVDATA.Opto.ActiveHitRate_LightOff = sum(LVDATA.Opto.HitTrialsArray(LightOffTrials(LightOffTrials<PosNaN),1)/nnz(~isnan(LVDATA.Opto.HitTrialsArray(LightOffTrials,1))));
LVDATA.Opto.ActiveHitRate_LightOn = sum(LVDATA.Opto.HitTrialsArray(LightOnTrials(LightOnTrials<PosNaN),1)/nnz(~isnan(LVDATA.Opto.HitTrialsArray(LightOnTrials,1))));
LVDATA.Opto.ActiveHitRate_LightAfter = sum(LVDATA.Opto.HitTrialsArray(AfterLightTrials(AfterLightTrials<PosNaN),1)/nnz(~isnan(LVDATA.Opto.HitTrialsArray(AfterLightTrials,1))));

%%% Binning light-stimulus signal to check in which position was it exactly
LVDATA.Opto.OptoStimOnMatrix = NaN(TRNu,LVDATA.BinNu);
% calculate the number of licks during each bin (lick/cm)
for i = 1:1:TRNu  % rows are trials
    for j = 1:1:LVDATA.BinNu  % columns (distance bins)  
        BinnedLight = mean(LVDATA.Opto.LightOnSignalDS(LVDATA.EnterIntoBinSampleInd(i,j):LVDATA.LeaveBinSampleInd(i,j)));
        if BinnedLight > 0.01 
            LVDATA.Opto.OptoStimOnMatrix(i,j) = 1;  % if there is a lick event, increase lick rate with one
        else
            LVDATA.Opto.OptoStimOnMatrix(i,j) = 0;
        end
    end
end
%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:TRNu,LVDATA.Opto.OptoStimOnMatrix) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
colormap(parula);
xlabel('Position on wheel (cm)');
ax = gca; ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(strcat(LVDATA.FileID,'-optical-stim'));
FileName = strcat('OptoStimHeatBin-',LVDATA.FileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

StimMax = max(mean(LVDATA.Opto.OptoStimOnMatrix));
LVDATA.Opto.optoStimStart = find(mean(LVDATA.Opto.OptoStimOnMatrix)>= StimMax/2,1)*LVDATA.BinSize;
LVDATA.Opto.optoStimEnd = find(mean(LVDATA.Opto.OptoStimOnMatrix)>= StimMax/2,1,'last')*LVDATA.BinSize;

end