function plotmeanROIActAllSessions(sData) % sData.imdata.MaoPC.LightOn or Off
% data: matrix in which rows are ROIs (mean activity during session in trials), columns are bins, binned Ca-activity dF/F values in cells (normalized to max activity)
% used normalized data (max = 1)% or do normalization
% GaussFilter = 5;


AllPlaceCellsPre = NaN(sData.imdata.nROIs,1);
for i=1:1:sData.imdata.nROIs %%% SET 
    if sum(any(sData1.imdata.MaoPC.PlaceCells(:,1)==i)) || sum(any(sData1.imdata.MaoPC.LightOn.PlaceCells(:,1)==i)) || sum(any(sData1.imdata.MaoPC.LightAfter.PlaceCells(:,1)==i))|| sum(any(sData2.imdata.MaoPC.PlaceCells(:,1)==i)|| sum(any(sData.imdata.MaoPC.PlaceCells(:,1)==i)))
        AllPlaceCellsPre(i,1) = i; 
    end
end
AllPlaceCells = AllPlaceCellsPre(~isnan(AllPlaceCellsPre));
ROINu= length(AllPlaceCells);
optoStimStart = sData1.behavior.opto.optoStimStart;
optoStimEnd = sData1.behavior.opto.optoStimEnd;

NumberOfDataSets = 5; %%% fill in
for i = 1:1:NumberOfDataSets
    Datasets{i} = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins); 
    DatasetName{i} = '';
end
Datasets{1,1} = sData1.imdata.binned.ROIsMeanAct_LightOffTrials; %SET
DatasetName{1} = 'S1-LaserOff';
Datasets{1,2} = sData1.imdata.binned.ROIsMeanAct_LightOnTrials; %SET
DatasetName{2} = 'S1-LaserOn';
Datasets{1,3} = sData1.imdata.binned.ROIsMeanAct_LightAfterTrials; %SET
DatasetName{3} = 'S1-LaserAfter';
Datasets{1,4} = sData2.imdata.binned.MeanRoiAct; %SET
DatasetName{4} = 'S2-LaserOn';
Datasets{1,5} = sData.imdata.binned.MeanRoiAct; %SET
DatasetName{5} = 'S3-LaserOff';


mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell\MeanActivitySortedAllSession');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell\MeanActivitySortedAllSession');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;



%SORTING based on peak of all sessions
PlaceDataSetMax(1:ROINu,1) =  max(Datasets{1,1}(AllPlaceCells,:),[],2);
PlaceDataSetMax(1:ROINu,2) =  max(Datasets{1,2}(AllPlaceCells,:),[],2);
PlaceDataSetMax(1:ROINu,3) =  max(Datasets{1,3}(AllPlaceCells,:),[],2);
PlaceDataSetMax(1:ROINu,4) =  max(Datasets{1,4}(AllPlaceCells,:),[],2);
PeakAmongSessions(1:ROINu,1) = AllPlaceCells; % place cells ROIS number
PeakAmongSessions(1:ROINu,2) = max(PlaceDataSetMax(),[],2); % max peak among all session 
for i = 1:1:ROINu
    PeakAmongSessions(i,3) = find(PlaceDataSetMax(i,:) == PeakAmongSessions(i,2)); % in which dataset (session) was the max
    PeakAmongSessions(i,4) = find(Datasets{1,PeakAmongSessions(i,3)}(AllPlaceCells(i),:)==PeakAmongSessions(i,2)); % in which bin was the max
end
SortingOrder = sortrows(PeakAmongSessions, 4); % first column is the sorted ROI order. plot in these order

for i = 1:1:NumberOfDataSets
    SortedDataPre = Datasets{1,i}(SortingOrder(:,1),:); % new matrix containing sorted data.
    Max = max(SortedDataPre,[],2);
    NormSortedData = SortedDataPre./Max;
    %PLOT FIGURE
    figure('Color','white'); 
    imagesc(Xaxis,1:ROINu,NormSortedData) %(1:number of bins;1:number of trials)
    c = colorbar; colormap(jet); caxis([0 1]);
    c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat('ROI activity sorted based on peak activity - ',DatasetName{i})); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-dff-SortedMeanAct-',DatasetName{i}); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end
    
end