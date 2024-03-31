function LVDATA = optoStimCalc(LVDATA,sData)  % OptoStimOnMatrix,optoStimStart,optoStimEnd

savePath = strcat(sData.sessionInfo.savePath,'\Behavior');

TRNu = LVDATA.TRNu;
PhotoStimLimit = (max(sData.daqdata.optoSignal) - min(sData.daqdata.optoSignal))/2; % size of photostim
LightOnSignal = zeros(size(sData.daqdata.optoSignal,1),1);
LightOnSignal(sData.daqdata.optoSignal > PhotoStimLimit) = 1;
LVDATA.Opto.LightOnSignalDS = LightOnSignal(LVDATA.FrameStartIndex);
%%% Binning light-stimulus signal to check in which position was it exactly
LVDATA.Opto.OptoStimOnMatrix = NaN(TRNu,LVDATA.BinNu);
% calculate the number of licks during each bin (lick/cm)
BinnedLightMatrix = NaN(TRNu,LVDATA.BinNu);
for i = 1:1:TRNu  % rows are trials
    for j = 1:1:LVDATA.BinNu  % columns (distance bins)  
        BinnedLight = mean(LVDATA.Opto.LightOnSignalDS(LVDATA.EnterIntoBinSampleInd(i,j):LVDATA.LeaveBinSampleInd(i,j)));
        BinnedLightMatrix(i,j) = BinnedLight;
        if BinnedLight > 0.5 
            LVDATA.Opto.OptoStimOnMatrix(i,j) = 1;  
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