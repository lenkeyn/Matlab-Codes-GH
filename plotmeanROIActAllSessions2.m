function plotmeanROIActAllSessions2(sData1,IsOptoSession,type) % sData.imdata.MaoPC.LightOn or Off
% data: matrix in which rows are ROIs (mean activity during session in trials), columns are bins, binned Ca-activity dF/F values in cells (normalized to max activity)
% used normalized data (max = 1)% or do normalization
% GaussFilter = 5; Type = 0;


AllPlaceCellsPre = NaN(sData1.imdata.nROIs,1);
for i=1:1:sData1.imdata.nROIs %%% SET!!!!! 
    if sum(any(sData1.imdata.MaoPC_dff.PlaceCells(:,1)==i)) || sum(any(sData2.imdata.MaoPC_dff.PlaceCells(:,1)==i)) || sum(any(sData3.imdata.MaoPC_dff.PlaceCells(:,1)==i)) 
        AllPlaceCellsPre(i,1) = i; 
    end
end
AllPlaceCells = AllPlaceCellsPre(~isnan(AllPlaceCellsPre));
ROINu= length(AllPlaceCells);
if IsOptoSession == 1
    optoStimStart = sData2.behavior.opto.optoStimStart;
    optoStimEnd = sData2.behavior.opto.optoStimEnd;
end

NumberOfDataSets = 3; %%% fill in!!!!
for i = 1:1:NumberOfDataSets
    Datasets{i} = NaN(sData1.imdata.nROIs,sData1.behavior.meta.nBins); 
    DatasetName{i} = '';
    Cues{i} = NaN(1,1);
    NormalizedSortedData{i} = NaN(1,1);
end
% fill up datasets DFF
if type == 0
    Datasets{1,1} = sData1.imdata.binned.MeanRoiAct; %SET
    Cues{1,1} = sData1.imdata.cues;
    DatasetName{1} = 'S1';
    Datasets{1,2} = sData2.imdata.binned.MeanRoiAct; %SET
    Cues{1,2} = sData2.imdata.cues;
    DatasetName{2} = 'S2';
    Datasets{1,3} = sData3.imdata.binned.MeanRoiAct; %SET
    DatasetName{3} = 'S3';
    Cues{1,3} = sData3.imdata.cues;
    %Datasets{1,4} = sData2.imdata.binned.MeanRoiAct; %SET
    %DatasetName{4} = 'S2-LaserOn';
    %Datasets{1,5} = sData.imdata.binned.MeanRoiAct; %SET
    %DatasetName{5} = 'S3-LaserOff';

    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'ComparingSessions\MeanActivitySortedAllSession-dff');
    savePath = strcat(sData1.sessionInfo.savePath,'\Imaging\ComparingSessions\MeanActivitySortedAllSession-dff');

elseif type == 1
% fill up datasets deconv
    Datasets{1,1} = sData1.imdata.binned.MeanRoiAct_Deconv; %SET
    Cues{1,1} = sData1.imdata.cues;
    DatasetName{1} = 'S1';
    Datasets{1,2} = sData2.imdata.binned.MeanRoiAct_Deconv; %SET
    Cues{1,2} = sData2.imdata.cues;
    DatasetName{2} = 'S2';
    Datasets{1,3} = sData3.imdata.binned.MeanRoiAct_Deconv; %SET
    DatasetName{3} = 'S3';

mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PlaceCell\MeanActivitySortedAllSession-deconv');
savePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PlaceCell\MeanActivitySortedAllSession-deconv');

elseif type == 2
% fill up datasets deconv
    Datasets{1,1} = sData1.imdata.binned.MeanRoiAct_FR; %SET
    DatasetName{1} = 'S1';
    Datasets{1,2} = sData2.imdata.binned.MeanRoiAct_FR; %SET
    DatasetName{2} = 'S2';
    Datasets{1,3} = sData3.imdata.binned.MeanRoiAct_FR; %SET
    DatasetName{3} = 'S3';

    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PlaceCell\MeanActivitySortedAllSession-FR');
    savePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PlaceCell\MeanActivitySortedAllSession-FR');
end

%SORTING based on peak of all sessions
PlaceDataSetMax(1:ROINu,1) =  max(Datasets{1,1}(AllPlaceCells,:),[],2);
PlaceDataSetMax(1:ROINu,2) =  max(Datasets{1,2}(AllPlaceCells,:),[],2);
PlaceDataSetMax(1:ROINu,3) =  max(Datasets{1,3}(AllPlaceCells,:),[],2); 
%PlaceDataSetMax(1:ROINu,4) =  max(Datasets{1,4}(AllPlaceCells,:),[],2);
PeakAmongSessions(1:ROINu,1) = AllPlaceCells; % place cells ROIS number
PeakAmongSessions(1:ROINu,2) = max(PlaceDataSetMax(),[],2); % max peak among all session 
for i = 1:1:ROINu
    PeakAmongSessions(i,3) = min(find(PlaceDataSetMax(i,:) == PeakAmongSessions(i,2))); % in which dataset (session) was the max
    PeakAmongSessions(i,4) = min(find(Datasets{1,PeakAmongSessions(i,3)}(AllPlaceCells(i),:)==PeakAmongSessions(i,2))); % in which bin was the max
end
SortingOrder = sortrows(PeakAmongSessions, 4); % first column is the sorted ROI order. plot in these order

Xaxis = sData1.behavior.meta.binSize:sData1.behavior.meta.binSize:sData1.behavior.meta.binSize*sData1.behavior.meta.nBins;
if type == 0
    fn = 'dff';
elseif type == 1
    fn = 'deconvolved';
    %%% smoothing for visualization
    for i = 1:1:NumberOfDataSets
        DeconvGauss = NaN(size(Datasets{1,i})); 
        for j = 1:1:length(DeconvGauss)
            SignalTemp1 = Datasets{1,i}(j,:);
            DeconvGauss(j,:) = smoothdata(SignalTemp1,'gaussian',5);
        end
        Datasets{1,i} = DeconvGauss;
    end
elseif type == 2
    fn = 'FR';
end

for i = 1:1:NumberOfDataSets
    SortedDataPre = Datasets{1,i}(SortingOrder(:,1),:); % new matrix containing sorted data.
    Max = max(SortedDataPre,[],2);
    NormSortedData = SortedDataPre./Max;
    NormalizedSortedData{i} = NormSortedData;
    %PLOT FIGURE
    figure('Color','white'); 
    imagesc(Xaxis,1:ROINu,NormSortedData) %(1:number of bins;1:number of trials)
    c = colorbar; colormap(jet); caxis([0 1]);
    c.Label.String = fn; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    if IsOptoSession == 1
        line([optoStimStart optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
        line([optoStimEnd optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    else
        line([Cues{1,i}.C1A Cues{1,i}.C1A],[0 length(AllPlaceCells)],'Color','black','LineStyle','--','LineWidth',2); hold on; line([Cues{1,i}.C1B Cues{1,i}.C1B],[0 length(AllPlaceCells)],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([Cues{1,i}.C2A Cues{1,i}.C2A],[0 length(AllPlaceCells)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([Cues{1,i}.C2B Cues{1,i}.C2B],[0 length(AllPlaceCells)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([Cues{1,i}.C3A Cues{1,i}.C3A],[0 length(AllPlaceCells)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([Cues{1,i}.C3B Cues{1,i}.C3B],[0 length(AllPlaceCells)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([Cues{1,i}.C4A Cues{1,i}.C4A],[0 length(AllPlaceCells)],'Color','black','LineStyle','--','LineWidth',2); hold on; line([Cues{1,i}.C4B Cues{1,i}.C4B],[0 length(AllPlaceCells)],'Color','black','LineStyle','--','LineWidth',2); hold on;
    end
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat('ROI activity sorted based on peak activity - ',DatasetName{i})); 
    FileName = strcat(sData1.sessionInfo.fileID,'-SortedMeanAct-',DatasetName{i}); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end
    
% Save sorting order
save(fullfile(savePath,strcat('AllSession-PC-sortingOrder',fn)),'SortingOrder');

end