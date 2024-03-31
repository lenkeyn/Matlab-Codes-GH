function LVDATA = OptoTrialSorting2(LVDATA,sData)
% TDMSDATA has been changed to sData.daqdata 20190221

LVDATA.Opto = struct;
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
    TrialLight = LVDATA.Opto.LightOnSignalDS(LVDATA.EnterIntoBinSampleInd(i,1):LVDATA.EnterIntoBinSampleInd(i+1,1));
    if sum(TrialLight)< numel(TrialLight)/100 % If there was light in the 1/100 the trial it should be considered as light-on trial
        LVDATA.Opto.LightOffTrials(i) = 1;
    else
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

%Lick
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

%%% Binning light-stimulus signal to check in which position was it exactly
LVDATA.Opto.OptoStimOnMatrix = NaN(TRNu,LVDATA.BinNu);
% calculate the number of licks during each bin (lick/cm)
for i = 1:1:TRNu  % rows are trials
    for j = 1:1:LVDATA.BinNu  % columns (distance bins)  
        BinnedLight = mean(LVDATA.Opto.LightOnSignalDS(LVDATA.EnterIntoBinSampleInd(i,j):LVDATA.LeaveBinSampleInd(i,j)));
        if BinnedLight > 0.5 
            LVDATA.Opto.OptoStimOnMatrix(i,j) = 1;  % if there is a lick event, increase lick rate with one
        else
            LVDATA.Opto.OptoStimOnMatrix(i,j) = 0;
        end
    end
end
%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:TRNu,LVDATA.Opto.OptoStimOnMatrix) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
colormap(jet);
xlabel('Position on wheel (cm)');
ax = gca; ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(strcat(LVDATA.FileID,'-optical-stim'));
FileName = strcat('OptoStimHeatBin-',LVDATA.FileID);
savefig(fullfile('C:\MATLAB\SAVE',FileName));
saveas(gcf,(fullfile('C:\MATLAB\SAVE',[FileName '.jpg'])));



end