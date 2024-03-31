function plotSortROIsMaxActHeatsDataPlaceOpto(sData,GaussFilter,type) % sData.imdata.MaoPC.LightOn or Off
% data: matrix in which rows are ROIs (mean activity during session in trials), columns are bins, binned Ca-activity dF/F values in cells (normalized to max activity)
% used normalized data (max = 1)% or do normalization
if type == 0
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell-dff\MeanActivitySorted');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell-dff\MeanActivitySorted');
elseif type == 1
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell-deconv\MeanActivitySorted');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell-deconv\MeanActivitySorted');    
end


%GaussFilter = 10; % cm 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;

if type == 0
    %light-off dff
    HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOffTrials(sData.imdata.MaoPC_dff.LightOff.PlaceCells,:);
    % sData.imdata.MaoPC.LightOff.PosTuningOrig(sData.imdata.MaoPC.LightOff.PlaceCells,:)
    % sData.imdata.binned.ROIsMeanAct_LightOffTrials(sData.imdata.MaoPC.LightOff.PlaceCells,:)
    ROINu= sData.imdata.MaoPC_dff.LightOff.placeROINu;
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
    c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'ROI activity sorted based on peak activity:Laser-off')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-dff-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end

%light-off deconv smoothed
if type ~= 0
    ROINu= sData.imdata.MaoPC_deconv.LightOff.placeROINu;
    GaussianSmoothed = NaN(length(sData.imdata.binned.LightOffTrials),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = sData.imdata.MaoPC_deconv.LightOff.PlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC_deconv.LightOff.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;

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
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'ROI activity sorted based on peak activity:Laser-off')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-deconvG-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end        

if type == 0 %light-on dff
    HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOnTrials(sData.imdata.MaoPC_dff.LightOn.PlaceCells,:);
    ROINu= sData.imdata.MaoPC_dff.LightOn.placeROINu;
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
    c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'ROI activity sorted based on peak activity:Laser-on')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-dff-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end

%light-on deconv
if type ~= 0
    ROINu= sData.imdata.MaoPC_deconv.LightOn.placeROINu;
    GaussianSmoothed = NaN(length(sData.imdata.binned.LightOnTrials),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = sData.imdata.MaoPC_deconv.LightOn.PlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC_deconv.LightOn.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;
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
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'ROI activity sorted based on peak activity:Laser-on')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-deconv SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end
  
if type == 0
    %light-after dff
    HeatMapData = sData.imdata.binned.ROIsMeanAct_LightAfterTrials(sData.imdata.MaoPC_dff.LightAfter.PlaceCells,:);
    ROINu= sData.imdata.MaoPC_dff.LightAfter.placeROINu;
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
    c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'ROI activity sorted based on peak activity:Laser-after')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserAfter-dFF-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end

%light-after deconv
if type ~= 0
    ROINu= sData.imdata.MaoPC_deconv.LightAfter.placeROINu;
    GaussianSmoothed = NaN(length(sData.imdata.binned.AfterLightTrials),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = sData.imdata.MaoPC_deconv.LightAfter.PlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC_deconv.LightAfter.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;
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
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title(strcat(sData.sessionInfo.fileID,'-ROI activity sorted based on peak activity:Laser-after')); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserAfter-deconv-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end

end