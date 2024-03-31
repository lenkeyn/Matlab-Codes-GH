function sData  = CaDataOptoSortingdFFsDataMoreProts(sData,FigVisible) 

%%% Set parameters:
datatype = 0; % datatype can be 0 = dFF or 1 = Deconvolved, 2= cumSpiking dataset
%Protocols = string({sData.stimProtocols.protocol});
Protocols = string({'none','stim 3 mW'}); %'stim 10 mW','stim 6mW','stim 3 mW ','stim 1 mW','stim 14-84 cm','stim 86-156 cm'

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActdFF-OptoMoreProt');
SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt'),'ROIActdFF');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt'),'MeanROIActdFF');
BinSize = sData.behavior.meta.binSize;
nBins = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging-1;
nROIs = sData.imdata.nROIs;
%nOptoProt = max(sData.behavior.optoMoreProts.OptoStimProtTrials); 
nOptoProt = length(unique(sData.behavior.optoMoreProts.OptoStimProtTrials)); % number of protocols


sData.imdata.binned.optoMoreProts = struct;
%Generating struct for storing data
for i = 1:1:nOptoProt
    for roi = 1:1:nROIs
        sData.imdata.binned.optoMoreProts.ROIsActdFF_Protocol{i,roi} = NaN(nTrials,nBins); % each Xth ROI's activity is the Xth a matrix in this cell array
        sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{i} = NaN(nROIs,nBins);
    end
end

%%% preparing data for plots
Xaxis = BinSize:BinSize:BinSize*nBins;
colors = jet(64);
colors(1,1:3) = 0.5; % set black when there is no signal


%{
%%% calculate opto-stim
optoStimMax = max(mean(sData.behavior.opto.OptoStimOnMatrix),[],2); % OptoStimOnMatrix: binarized matrix (0-1) if opto stimulus is on or off 
meanOptoStim = zeros(1,sData.behavior.meta.nBins);
meanOptoStim((mean(sData.behavior.opto.OptoStimOnMatrix)) > optoStimMax/4) = 1;
optoStimStart = min(find(diff(meanOptoStim) == 1)*sData.behavior.meta.binSize);
optoStimEnd = max(find(diff(meanOptoStim) == -1)*sData.behavior.meta.binSize);
if  isempty(optoStimEnd) || optoStimEnd < optoStimStart
    optoStimEnd = 157;
end
sData.behavior.opto.optoStimStart = optoStimStart;
sData.behavior.opto.optoStimEnd = optoStimEnd;
%}

optoStimStart = sData.behavior.opto.optoStimStart;
optoStimEnd = sData.behavior.opto.optoStimEnd;

for roi = 1:1:sData.imdata.nROIs
    
    if datatype == 0
        data = sData.imdata.binned.RoidFF{1,roi};
    elseif datatype == 1
        data = sData.imdata.binned.RoiDeconvolved{1,roi};
    elseif datatype == 2
        data = sData.imdata.binned.RoicumSpikes{1,roi};    
    end
    
    
    for i = 1:1:nOptoProt 
        k = find(sData.behavior.optoMoreProts.OptoStimProtTrials == i);
        sData.imdata.binned.optoMoreProts.ROIsActdFF_Protocol{i,roi}(k,:) = data(k,:);
        sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{i}(roi,:) = nanmean(data(k,:),1);
    end

    Cmax = max(max(data));
    Cmin = min(min(data));
    if datatype == 0
        Min = Cmin;
        text = 'dFF';
    elseif datatype == 1
        Min = -Cmax/10;  
        text = 'Deconvolved signal';
    elseif datatype == 2
        Min = -Cmax/10;  
        text = 'Spikerate';
    end 
     
    
   %%% PLOT
    
    figure('Color','white','Position',[100 100 900 600],'visible',FigVisible); 
    Ymax = 0.0001;
    % heatplots
    for j=1:1:nOptoProt
        subplot(3,nOptoProt,[j, j+nOptoProt]); 
        imagesc(Xaxis,1:nTrials,sData.imdata.binned.optoMoreProts.ROIsActdFF_Protocol{j,roi}); 
        if j>1
            line([optoStimStart optoStimStart],[0 nTrials+1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials+1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        end
        colormap(colors); caxis([Min Cmax]); %Cmin/10 for dff  % -Cmax/10 for spikerate
        xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials','FontSize',12);
        title(Protocols(j),'FontSize',12);
        if Ymax < max(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(roi,:)) % looking for Ymax value for the next plot
            Ymax = max(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(roi,:));
        end
    end
    % mean act plots
    for j=1:1:nOptoProt
        subplot(3,nOptoProt,j+2*nOptoProt); 
        plot(Xaxis,sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(roi,:)); 
        if j>1
            line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        end
        axis([0 160 0 Ymax]); % ceil(Ymax)
        xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning','FontSize',12);
    end 
    suptitleN(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi),'-dFF'));    
    FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-optoMoreProt'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt\ROIActdFF'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt\ROIActdFF'),[FileName '.jpg'])));
    
    
    %%% PLOT mean ROI activity in different protocols in one plot
    figure('Color','white','visible',FigVisible)
    for j=1:1:nOptoProt
        plot(Xaxis,smoothdata(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(roi,:),'gaussian',5),'LineWidth',1.5); hold on;
    end
    line([optoStimStart optoStimStart],[-Ymax/5 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[-Ymax/5 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    %axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Mean position tuning','FontSize',12);
    legend(Protocols,'Location','south','FontSize',12); 
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi),'-dFF'));
    FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-optoMoreProtMean-dFF-G5'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt\MeanROIActdFF'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiActdFF-OptoMoreProt\MeanROIActdFF'),[FileName '.jpg'])));
    
    % close window if there are too many 
    if rem(roi,25) == 0 
        close all;
    end

end
close all;

%%% PLOT mean ROI activity in different protocols in one plot
figure('Color','white')
Ymax = 0;
for j=1:1:nOptoProt
    plot(Xaxis,nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(:,:),1),'LineWidth',1.5); hold on; %smoothdata(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(:,:),1)
    if max(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(:,:),1),[],2) > Ymax
        Ymax = max(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanActdFF_Protocol{1,j}(:,:),1),[],2);
    end
end
line([optoStimStart optoStimStart],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
%axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Mean position tuning','FontSize',12);
legend(Protocols,'Location','South','FontSize',12); 
title(strcat(sData.sessionInfo.fileID,'-all-ROIs'));
FileName = strcat(sData.sessionInfo.fileID,'-optoMoreProtMeanAllROIsdFF'); 
savefig(fullfile(SavePath,FileName));
saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
 
end