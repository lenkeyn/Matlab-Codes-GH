function sData = plotSortROIsMaxActHeatsData(sData) % PC_Mao1.SpikeRate_AllTrialsNorm
% parameters: ;
% data: matrix in which rows are ROIs (mean activity during session in trials), columns are bins, binned Ca-activity dF/F values in cells (normalized to max activity)
% used normalized data (max = 1)
% or do normalization:
[ROINu, ~] = size(sData.imdata.binned.ROIsMeanAct_LightOnTrials);
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;

%light-off
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOffTrials;
Max = max(HeatMapData,[],2);
NormHeatMapData = HeatMapData./Max;
MaxData = max(NormHeatMapData, [], 2);
MaxDataBin = NaN(ROINu,2);
MaxDataBin(:,2) = (1:ROINu);
for i = 1:1:ROINu
    MaxDataBin(i,1) = find(NormHeatMapData(i,:) == MaxData(i));
end
SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
SortedData1 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,1:ROINu,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'Normalized Spike Rate'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('Place-cell activity sorted based on peak activity - Light-off'); 
FileName = strcat(sData.sessionInfo.fileID,'-LightOff-SortedMeanAct'); 
savefig(fullfile(sData.sessionInfo.savePath,FileName));
saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));
        

%light-on
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOnTrials;
Max = max(HeatMapData,[],2);
NormHeatMapData = HeatMapData./Max;
MaxData = max(NormHeatMapData, [], 2);
MaxDataBin = NaN(ROINu,2);
MaxDataBin(:,2) = (1:ROINu);
for i = 1:1:ROINu
    MaxDataBin(i,1) = find(NormHeatMapData(i,:) == MaxData(i));
end
SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
SortedData2 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,1:ROINu,SortedData2) %(1:number of bins;1:number of trials) 
c = colorbar; colormap(jet); caxis([0 1]); hold on;
c.Label.String = 'Normalized Spike Rate'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('Place-cell activity sorted based on peak activity - Light-on'); 
FileName = strcat(sData.sessionInfo.fileID,'-LightOn-SortedMeanAct'); 
savefig(fullfile(sData.sessionInfo.savePath,FileName));
saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));
        

%light-after
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightAfterTrials;
Max = max(HeatMapData,[],2);
NormHeatMapData = HeatMapData./Max;
MaxData = max(NormHeatMapData, [], 2);
MaxDataBin = NaN(ROINu,2);
MaxDataBin(:,2) = (1:ROINu);
for i = 1:1:ROINu
    MaxDataBin(i,1) = find(NormHeatMapData(i,:) == MaxData(i));
end
SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
SortedData3 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,1:ROINu,SortedData3) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'Normalized Spike Rate'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('Place-cell activity sorted based on peak activity - After-Light'); 
FileName = strcat(sData.sessionInfo.fileID,'-LightAfter-SortedMeanAct'); 
savefig(fullfile(sData.sessionInfo.savePath,FileName));
saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));

% calculation of difference of light-off and light-on trials
% SortedData1 is the light-off data of mean roi act, SortedData4 will be the light-on trials data, normalized to light-off trials
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOnTrials;
NormHeatMapData = HeatMapData./Max;
MaxData = max(NormHeatMapData, [], 2);
MaxDataBin = NaN(ROINu,2);
MaxDataBin(:,2) = (1:ROINu);
for i = 1:1:ROINu
    MaxDataBin(i,1) = find(NormHeatMapData(i,:) == MaxData(i));
end
SortedData4 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
sData.imdata.binned.ROIsOptoEffectHeat = SortedData1-SortedData4;
figure('Color','white'); 
imagesc(Xaxis,1:ROINu,sData.imdata.binned.ROIsOptoEffectHeat)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'Normalized Spike Rate'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('Optical stimulation effect on activity (light-off vs light-on trials)'); 
FileName = strcat(sData.sessionInfo.fileID,'OptoEffectROIHeat-LightOn-Off'); 
savefig(fullfile(sData.sessionInfo.savePath,FileName));
saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));        

end