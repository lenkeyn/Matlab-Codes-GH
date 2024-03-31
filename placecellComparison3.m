function sData = placecellComparison3(sData,type,GaussFilter)
%type = 1; %0 use dff data, 1 use deconvolved data, 2 use firing rate data
% GaussFilter = 5;

AllPlaceCellsPre = NaN(sData.imdata.nROIs,1);
for i=1:1:sData.imdata.nROIs 
     %{
    if sum(any(sData.imdata.MaoOpto_dff{1,1}.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoOpto_dff{1,2}.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoOpto_dff{1,3}.PlaceCells(:,1)==i))% || sum(any(sData.imdata.MaoOpto_dff{1,4}.PlaceCells(:,1)==i)) 
        AllPlaceCellsPre(i,1) = i;  
    end
    %}
    if type == 0 
        if sum(any(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceCells(:,1)==i))
            AllPlaceCellsPre(i,1) = i;  
        end
    end
    if type ~= 0 
        if sum(any(sData.imdata.MaoPC_Opto_deconv.OptoOff.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_Opto_deconv.OptoOn.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_Opto_deconv.OptoAfter.PlaceCells(:,1)==i))
            AllPlaceCellsPre(i,1) = i;  
        end
    end
    %}
    %{ 
    if sum(any(sData.imdata.MaoPC.OptoOff.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC.OptoOn.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC.OptoAfter.PlaceCells(:,1)==i))
        AllPlaceCellsPre(i,1) = i;  
    end
        %sData1.imdata.MaoPC_dff.PlaceCells
    
   %}
    %{
    if sum(any(sData.imdata.MaoPC_dff.OptoOff.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_dff.OptoOn.PlaceCells(:,1)==i)) || sum(any(sData.imdata.MaoPC_dff.OptoAfter.PlaceCells(:,1)==i)) %|| sum(any(sData.imdata.MaoOpto_dff{1,4}.PlaceCells(:,1)==i))  
        AllPlaceCellsPre(i,1) = i;  
    end
     %}
end
AllPlaceCells = AllPlaceCellsPre(~isnan(AllPlaceCellsPre));

if type == 0 
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MeanActivitySortedComparison-dff');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\MeanActivitySortedComparison-dff');
else
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MeanActivitySortedComparison-deconv');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\MeanActivitySortedComparison-deconv');
end

Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
ROINu= length(AllPlaceCells);
Yaxis = 1:1:ROINu; 


if type == 0
    %%% Opto off , Opto on, Opto after trials
    Max = NaN(length(AllPlaceCells),3);
    Max(:,1) = max(sData.imdata.binned.ROIsMeanAct_OptoOffTrials(AllPlaceCells,:),[],2);
    Max(:,2) = max(sData.imdata.binned.ROIsMeanAct_OptoOnTrials(AllPlaceCells,:),[],2);
    Max(:,3) = max(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(AllPlaceCells,:),[],2);
    MaxAct = max(Max,[],2);

    %laser-off dff
    % nomrlize data
    HeatMapData = sData.imdata.binned.ROIsMeanAct_OptoOffTrials(AllPlaceCells,:); % sData.imdata.MaoPC.OptoOff.PosTuningOrig(sData.imdata.MaoPC.OptoOff.PlaceCells,:)
    NormHeatMapData = HeatMapData./MaxAct;
    % sorting
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
end

% deconv smoothed calculation
if type ~= 0
      
    GaussianSmoothedOptoOff = NaN(length(AllPlaceCells),sData.behavior.meta.nBins); %length(sData.imdata.binned.OptoOffTrials)
    GaussianSmoothedOptoOn = NaN(length(AllPlaceCells),sData.behavior.meta.nBins);
    GaussianSmoothedOptoAfter = NaN(length(AllPlaceCells),sData.behavior.meta.nBins);
    for i = 1:1:ROINu
        j = AllPlaceCells(i);
        GaussianSmoothedOptoOff(i,:) = smoothdata(sData.imdata.MaoPC_Opto_deconv.OptoOff.PosTuningOrig(j,:),'gaussian',GaussFilter);
        GaussianSmoothedOptoOn(i,:) = smoothdata(sData.imdata.MaoPC_Opto_deconv.OptoOn.PosTuningOrig(j,:),'gaussian',GaussFilter);
        GaussianSmoothedOptoAfter(i,:) = smoothdata(sData.imdata.MaoPC_Opto_deconv.OptoAfter.PosTuningOrig(j,:),'gaussian',GaussFilter);
    end
    
    MaxDeconv = NaN(length(AllPlaceCells),3);
    MaxDeconv(:,1) = max(GaussianSmoothedOptoOff,[],2);
    MaxDeconv(:,2) = max(GaussianSmoothedOptoOn,[],2);
    MaxDeconv(:,3) = max(GaussianSmoothedOptoAfter,[],2);
    MaxDeconvAct = max(MaxDeconv,[],2);
    
    % plot laser-off deconv smoothed
    HeatMapData = GaussianSmoothedOptoOff;
    %Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./MaxDeconvAct;
    MaxData = max(NormHeatMapData, [], 2);
    MaxDataBin = NaN(ROINu,2);
    MaxDataBin(:,2) = (1:ROINu);
    for i = 1:1:ROINu
        MaxDataBin(i,1) = find(NormHeatMapData(i,:) == MaxData(i),1);
    end
    SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
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

if type == 0
    %Opto-on dff
    HeatMapData = sData.imdata.binned.ROIsMeanAct_OptoOnTrials(AllPlaceCells,:);
    %Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./MaxAct;
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
end

%Opto-on deconv
if type ~= 0
    HeatMapData = GaussianSmoothedOptoOn;
    %Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./MaxDeconvAct;
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
    
if type == 0
    %Opto-after dff
    HeatMapData = sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(AllPlaceCells,:);
    %Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./MaxAct;
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
end

%Opto-after deconv
if type ~= 0
    HeatMapData = GaussianSmoothedOptoAfter;
    %Max = max(HeatMapData,[],2);
    NormHeatMapData = HeatMapData./MaxDeconvAct;
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
end

% PLOT potition tuning of all place cells
if type == 0
    %Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    MeanActLaserOff = nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials(AllPlaceCells,:));
    MeanActLaserOn = nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials(AllPlaceCells,:));
    MeanActLaserAfter = nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(AllPlaceCells,:));
end
if type ~= 0
    MeanActLaserOff = nanmean(SortedData1,1);
    MeanActLaserOn = nanmean(SortedData2,1);
    MeanActLaserAfter = nanmean(SortedData3,1);
end

    Ymax = max([max(MeanActLaserOff) max(MeanActLaserOn) max(MeanActLaserAfter)])*1.05;
    Ymin = min(MeanActLaserOff)*0.9;

    figure('Color','white','Position',[200 200 1000 200]);
    % laser-off mean act
    subplot(1,3,1); 
    plot(Xaxis,MeanActLaserOff) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('Laser-off'); 

    % laser-on mean act
    subplot(1,3,2); 
    plot(Xaxis,MeanActLaserOn)  %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOnTrials(roi,1:BinNu)
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 Ymin Ymax]);
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('Laser-on');

    % laser-after mean act
    subplot(1,3,3); 
    plot(Xaxis,MeanActLaserAfter)   %Xaxis,CADATA.Opto.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu)
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 Ymin Ymax]);
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('Laser-after');
    FileName = strcat(sData.sessionInfo.fileID,'-optoAllPlaceCells'); 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    if type == 0
    sData.imdata.MaoPC_Sorting_dFF.AllPlaceCells = AllPlaceCells;
    sData.imdata.MaoPC_Sorting_dFF.SortingOrder = SortingOrder;
    sData.imdata.MaoPC_Sorting_dFF.AllPC.MeanActLaserOn = MeanActLaserOn;
    sData.imdata.MaoPC_Sorting_dFF.AllPC.MeanActLaserOn = MeanActLaserOff;
    sData.imdata.MaoPC_Sorting_dFF.AllPC.MeanActLaserAfter = MeanActLaserAfter;
    else
    sData.imdata.MaoPC_Sorting_deconv.AllPlaceCells = AllPlaceCells;
    sData.imdata.MaoPC_Sorting_deconv.SortingOrder = SortingOrder;
    sData.imdata.MaoPC_Sorting_deconv.AllPC.MeanActLaserOn = MeanActLaserOn;
    sData.imdata.MaoPC_Sorting_deconv.AllPC.MeanActLaserOn = MeanActLaserOff;
    sData.imdata.MaoPC_Sorting_deconv.AllPC.MeanActLaserAfter = MeanActLaserAfter;    
    end


    %mean activity of all ROIs Opto-off-on-after in one plot
    figure('Color','white','Position',[100 100 500 400]);
    plot(Xaxis,MeanActLaserOff,'LineStyle','-','LineWidth',1) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    hold on
    plot(Xaxis,MeanActLaserOn,'LineStyle','-','LineWidth',1)
    hold on
    plot(Xaxis,MeanActLaserAfter,'LineStyle','-','LineWidth',1)
    hold on
    line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 sData.behavior.wheelLapImaging],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    legend('Laser-off', 'Laser-on', 'After-laser','Location','south');
    title('Position tuning of all place cells');
    FileName = '-optoAllPlaceCells2'; 
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

 %{
%%% multiple protocols in one session
Max = NaN(length(AllPlaceCells),4);
Max(:,1) = max(sData.imdata.MaoOpto_dff{1, 1}.PosTuning(AllPlaceCells,:),[],2);
Max(:,2) = max(sData.imdata.MaoOpto_dff{1, 2}.PosTuning(AllPlaceCells,:),[],2);
Max(:,3) = max(sData.imdata.MaoOpto_dff{1, 3}.PosTuning(AllPlaceCells,:),[],2);
%Max(:,4) = max(sData.imdata.MaoOpto_dff{1, 4}.PosTuning(AllPlaceCells,:),[],2);
MaxAct = max(Max,[],2);


%laser-off control dff (protocol #1 )
% nomrlize data
HeatMapData = sData.imdata.MaoOpto_dff{1, 1}.PosTuning(AllPlaceCells,:); % sData.imdata.MaoPC.OptoOff.PosTuningOrig(sData.imdata.MaoPC.OptoOff.PlaceCells,:)
NormHeatMapData = HeatMapData./MaxAct;
% sorting
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
title('ROI activity sorted based on peak activity - Prot1'); 
FileName = strcat(sData.sessionInfo.fileID,'-Prot1-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%protocol #2 
HeatMapData =  sData.imdata.MaoOpto_dff{1, 2}.PosTuning(AllPlaceCells,:);
NormHeatMapData = HeatMapData./MaxAct;
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
title('ROI activity sorted based on peak activity - Prot2'); 
FileName = strcat(sData.sessionInfo.fileID,'-Prot2-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%protocol #3 
HeatMapData =  sData.imdata.MaoOpto_dff{1, 3}.PosTuning(AllPlaceCells,:);
NormHeatMapData = HeatMapData./MaxAct;
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
title('ROI activity sorted based on peak activity - Prot3'); 
FileName = strcat(sData.sessionInfo.fileID,'-Prot3-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%protocol #4 
HeatMapData =  sData.imdata.MaoOpto_dff{1, 4}.PosTuning(AllPlaceCells,:);
NormHeatMapData = HeatMapData./MaxAct;
SortedData4 = NormHeatMapData(SortingOrder(:,2),:); % new matrix containing sorted data.  % NormHeatMapData
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData4) %(1:number of bins;1:number of trials) 
c = colorbar; colormap(jet); caxis([0 1]); hold on;
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 ROINu],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('ROI activity sorted based on peak activity - Prot4'); 
FileName = strcat(sData.sessionInfo.fileID,'-Prot4-dff-SortedMeanAct'); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%mean activity of all ROIs in one plot
figure('Color','white','Position',[100 100 500 400]);
Ymax = max([mean(SortedData1) mean(SortedData2) mean(SortedData3)])*1.05;
plot(Xaxis,mean(SortedData1),'LineStyle','-','LineWidth',1) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
hold on
plot(Xaxis,mean(SortedData2),'LineStyle','-','LineWidth',1)
hold on
plot(Xaxis,mean(SortedData3),'LineStyle','-','LineWidth',1)
hold on
%plot(Xaxis,mean(SortedData4),'LineStyle','-','LineWidth',1)
%hold on
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning');
legend('Laser-off', 'Laser-on 1 mW', 'Laser-on 0.3 mW','Location','south');
%legend('Laser-off', 'Laser-on 14-84 cm', 'Laser-on 86-156 cm','Location','southwest'); %'Laser-on 10 mW', 'Laser-on 3 mW', 'Laser-on 1 mW',
%legend('Laser-off', 'Laser-on 14-84cm', 'Laser-on 86-156cm','Location','southwest');
title('Position tuning of all place cells');
FileName = '-optoAllPlaceCellsMoreProt'; 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%}

end

