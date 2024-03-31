function placecellComparison(sData,type,GaussFilter)
%type = 1; %0 use dff data, 1 use deconvolved data, 2 use firing rate data
% GaussianFilter = 5;

AllPlaceCellsPre = NaN(sData.imdata.nROIs,1);
for i=1:1:sData.imdata.nROIs
    if sum(any(sData.imdata.MaoPC.LightOff.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC.LightOn.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC.LightAfter.PlaceCells(:,1)==i))
        AllPlaceCellsPre(i,1) = i; 
    end
end
AllPlaceCells = AllPlaceCellsPre(~isnan(AllPlaceCellsPre));

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell\MeanActivitySortedComparison');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell\MeanActivitySortedComparison');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
ROINu= length(AllPlaceCells);
Yaxis = 1:1:ROINu; 

Max = NaN(length(AllPlaceCells),3);
Max(:,1) = max(sData.imdata.binned.ROIsMeanAct_LightOffTrials(AllPlaceCells,:),[],2);
Max(:,2) = max(sData.imdata.binned.ROIsMeanAct_LightOnTrials(AllPlaceCells,:),[],2);
Max(:,3) = max(sData.imdata.binned.ROIsMeanAct_LightAfterTrials(AllPlaceCells,:),[],2);
MaxAct = max(Max,[],2);

%laser-off dff
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOffTrials(AllPlaceCells,:); % sData.imdata.MaoPC.LightOff.PosTuningOrig(sData.imdata.MaoPC.LightOff.PlaceCells,:)
NormHeatMapData = HeatMapData./MaxAct;
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
imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('ROI activity sorted based on peak activity - Laser-off'); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%laser-off deconv smoothed
if type ~= 0
    GaussianSmoothed = NaN(length(sData.imdata.binned.LightOffTrials),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = AllPlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC.LightOff.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;
    Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./Max;
    SortedData1 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.
    %PLOT FIGURE
    figure('Color','white'); 
    imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
    c = colorbar; colormap(jet); caxis([0 1]);
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title('ROI activity sorted based on peak activity - Laser-off'); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-deconvG-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end        

%light-on dff
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightOnTrials(AllPlaceCells,:);
Max = max(HeatMapData,[],2);
NormHeatMapData = HeatMapData./Max;
SortedData2 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData2) %(1:number of bins;1:number of trials) 
c = colorbar; colormap(jet); caxis([0 1]); hold on;
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('ROI activity sorted based on peak activity - Laser-on'); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%light-on deconv
if type ~= 0
    GaussianSmoothed = NaN(length(AllPlaceCells),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = AllPlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC.LightOn.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;
    Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./Max;
    SortedData2 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
    %PLOT FIGURE
    figure('Color','white'); 
    imagesc(Xaxis,Yaxis,SortedData2) %(1:number of bins;1:number of trials) 
    c = colorbar; colormap(jet); caxis([0 1]); hold on;
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title('ROI activity sorted based on peak activity - Laser-on'); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-deconv SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end
    
%light-after dff
HeatMapData = sData.imdata.binned.ROIsMeanAct_LightAfterTrials(AllPlaceCells,:);
Max = max(HeatMapData,[],2);
NormHeatMapData = HeatMapData./Max;
SortedData3 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData3) %(1:number of bins;1:number of trials) 
c = colorbar; colormap(jet); caxis([0 1]); hold on;
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('ROI activity sorted based on peak activity - Laser-after'); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserAfter-dFF-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%light-after deconv
if type ~= 0
    GaussianSmoothed = NaN(length(AllPlaceCells),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = AllPlaceCells(i);
        GaussianSmoothed(i,:) = smoothdata(sData.imdata.MaoPC.LightAfter.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    HeatMapData = GaussianSmoothed;
    Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./Max;
    SortedData3 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
    %PLOT FIGURE
    figure('Color','white'); 
    imagesc(Xaxis,Yaxis,SortedData3) %(1:number of bins;1:number of trials) 
    c = colorbar; colormap(jet); caxis([0 1]); hold on;
    c.Label.String = 'Smoothed deconvolved data'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
    line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
    set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
    title('ROI activity sorted based on peak activity - Laser-after'); 
    FileName = strcat(sData.sessionInfo.fileID,'-LaserAfter-deconv-SortedMeanAct'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% PLOT potition tuning of all place cells    

%Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
MeanActLaserOff = nanmean(sData.imdata.binned.ROIsMeanAct_LightOffTrials(AllPlaceCells,:));
MeanActLaserOn = nanmean(sData.imdata.binned.ROIsMeanAct_LightOnTrials(AllPlaceCells,:));
MeanActLaserAfter = nanmean(sData.imdata.binned.ROIsMeanAct_LightAfterTrials(AllPlaceCells,:));
Ymax = max(MeanActLaserOff)*1.1;
Ymin = min(MeanActLaserOff)*0.9;

figure('Color','white','Position',[200 200 1000 200]);
% laser-off mean act
subplot(1,3,1); 
plot(Xaxis,MeanActLaserOff) %Xaxis,CADATA.Opto.ROIsMeanAct_LightOffTrials(roi,1:BinNu)
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 Ymin Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning');
title('Laser-off'); 

% laser-on mean act
subplot(1,3,2); 
plot(Xaxis,MeanActLaserOn)  %Xaxis,CADATA.Opto.ROIsMeanAct_LightOnTrials(roi,1:BinNu)
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 Ymin Ymax]);
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning');
title('Laser-on');

% laser-after mean act
subplot(1,3,3); 
plot(Xaxis,MeanActLaserAfter)   %Xaxis,CADATA.Opto.ROIsMeanAct_LightAfterTrials(roi,1:BinNu)
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 Ymin Ymax]);
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning');
title('Laser-after');
FileName = strcat(sData.sessionInfo.fileID,'-optoAllPlaceCells'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

sData.imdata.MaoPC.AllPlaceCells = AllPlaceCells;
sData.imdata.MaoPC.SortingOrder = SortingOrder;
sData.imdata.MaoPC.AllPC.MeanActLaserOn = MeanActLaserOn;
sData.imdata.MaoPC.AllPC.MeanActLaserOn = MeanActLaserOff;
sData.imdata.MaoPC.AllPC.MeanActLaserAfter = MeanActLaserAfter;

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end