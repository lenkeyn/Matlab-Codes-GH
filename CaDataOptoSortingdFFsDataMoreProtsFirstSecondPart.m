function sData  = CaDataOptoSortingdFFsDataMoreProtsFirstSecondPart(sData,type,FigVisible) % type can be 0 = RoidFF or 1 = RoidFF-SR, 2= RoidFF-SR-LP

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct-OptoMoreProt');
SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt'),'ROIAct');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt'),'MeanROIAct');
BinSize = sData.behavior.meta.binSize;
nBins = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging-1;
nROIs = sData.imdata.nROIs;
nOptoProt = max(sData.behavior.optoMoreProts.OptoStimProtTrials); 

%Protocols = string({sData.stimProtocols.protocol});
Protocols = string({'none','stim14-84','stim86-156'}); %'stim10','stim3','stim1','stim07'

sData.imdata.binned.optoMoreProts = struct;
%Generating struct for storing data
for i = 1:1:nOptoProt
    for roi = 1:1:nROIs
    sData.imdata.binned.optoMoreProts.ROIsAct_Protocol{i,roi} = NaN(nTrials,nBins); % each Xth ROI's activity is the Xth a matrix in this cell array
    sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{i} = NaN(nROIs,nBins);
    end
end

%%% preparing data for plots
Xaxis = BinSize:BinSize:BinSize*nBins;
colors = jet(64);
colors(1,1:3) = 0.5; % set black when there is no signal

%%% calculate opto-stim
optoStimMax = max(mean(sData.behavior.opto.OptoStimOnMatrix),[],2);
meanOptoStim = zeros(1,sData.behavior.meta.nBins);
meanOptoStim((mean(sData.behavior.opto.OptoStimOnMatrix)) > optoStimMax/4) = 1;
optoStimStart = find(diff(meanOptoStim) == 1)*sData.behavior.meta.binSize;
optoStimEndpre = find(diff(meanOptoStim) == -1)*sData.behavior.meta.binSize;
optoStimEnd(1) = optoStimEndpre(2);
optoStimEnd(2) = 157;
%sData.behavior.opto.optoStimStart = optoStimStart;
%sData.behavior.opto.optoStimEnd = optoStimEnd;

for roi = 1:1:sData.imdata.nROIs
    
    if type == 0
        data = sData.imdata.binned.RoidFF{1,roi};
    elseif type == 1
        data = sData.imdata.binned.RoiDeconvolved{1,roi};
    elseif type == 2
        data = sData.imdata.binned.RoiSpikeRate{1,roi};    
    end
    
    
    for i = 1:1:nOptoProt %
        k = find(sData.behavior.optoMoreProts.OptoStimProtTrials == i);
        sData.imdata.binned.optoMoreProts.ROIsAct_Protocol{i,roi}(k,:) = data(k,:);
        sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{i}(roi,:) = nanmean(data(k,:),1);
    end

    Cmax = max(max(data));
    Cmin = min(min(data));
    if type == 0
        Min = Cmin;
        text = 'dFF';
    elseif type == 1
        Min = -Cmax/10;  
        text = 'Spikerate';
    end 
     
   %%% PLOT
    
    figure('Color','white','Position',[100 100 1500 800],'visible',FigVisible); 
    Ymax = 0.000000001;
    % heatplots
    for j=1:1:nOptoProt
        subplot(3,nOptoProt,[j, j+nOptoProt]); 
        imagesc(Xaxis,1:nTrials,sData.imdata.binned.optoMoreProts.ROIsAct_Protocol{j,roi}); 
        if j>1
            line([optoStimStart(j-1) optoStimStart(j-1)],[0 nTrials+1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd(j-1) optoStimEnd(j-1)],[0 nTrials+1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        end
        colormap(colors); caxis([Min Cmax]); %Cmin/10 for dff  % -Cmax/10 for spikerate
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials');
        title(Protocols(j),'FontSize',7);
        if Ymax < max(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(roi,:)) % looking for Ymax value for the next plot
            Ymax = max(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(roi,:));
        end
    end
    % mean act plots
    for j=1:1:nOptoProt
        subplot(3,nOptoProt,j+2*nOptoProt); 
        plot(Xaxis,sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(roi,:)); 
        if j > 1
            line([optoStimStart(j-1) optoStimStart(j-1)],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd(j-1) optoStimEnd(j-1)],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        end
        axis([0 160 0 Ymax]); % ceil(Ymax)
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning');
        
    end 
    suptitleN(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi)));    
    FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-optoMoreProt'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt\ROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt\ROIAct'),[FileName '.jpg'])));
    
    
    %%% PLOT mean ROI activity in different protocols in one plot
    figure('Color','white','visible',FigVisible)
    for j=1:1:nOptoProt
        plot(Xaxis,smoothdata(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(roi,:),'gaussian',5),'LineWidth',1.5); hold on;
        if j>1
            line([optoStimStart(j-1) optoStimStart(j-1)],[-Ymax/5 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd(j-1) optoStimEnd(j-1)],[-Ymax/5 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        end
    end
    %axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Mean position tuning');
    legend(Protocols,'Location','southoutside'); 
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi)));
    FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-optoMoreProtMean'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt\MeanROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-OptoMoreProt\MeanROIAct'),[FileName '.jpg'])));
    
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
    plot(Xaxis,nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(:,:),1),'LineWidth',1.5); hold on; %smoothdata(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(:,:),1)
    if max(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(:,:),1),[],2) > Ymax
        Ymax = max(nanmean(sData.imdata.binned.optoMoreProts.ROIsMeanAct_Protocol{1,j}(:,:),1),[],2);
    end
end
line([optoStimStart(1) optoStimStart(1)],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd(1) optoStimEnd(1)],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
line([optoStimStart(2) optoStimStart(2)],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd(2) optoStimEnd(2)],[0 Ymax*1.1],'Color','red','LineStyle','-','LineWidth',1); hold on; 
%axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Mean position tuning');
legend(Protocols,'Location','South'); 
title(strcat(sData.sessionInfo.fileID,'-all-ROIs'));
FileName = strcat(sData.sessionInfo.fileID,'-optoMoreProtMean'); 
savefig(fullfile(SavePath,FileName));
saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
 
end