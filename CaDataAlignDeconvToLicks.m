function sData  = CaDataAlignDeconvToLicks(sData,FigVisible) 

% FigVisible = 'on';
%params:
DataShownAfterStimMs = 1000; % data shown after stimulus
DataShownBeforeStimMs = 1000; % data shown before stimulus 
DataShownAfterStimSamples = round(DataShownAfterStimMs/sData.imdata.meta.fps);
DataShownBeforeStimSamples = round(DataShownBeforeStimMs/sData.imdata.meta.fps);
DataShownLength = DataShownAfterStimSamples+DataShownBeforeStimSamples+1;

% set
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'AlignDffToLicks');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks'),'ROIAct');
%SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim');
%nTrials = sData.behavior.wheelLapImaging-1;
%nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;

%create struct
sData.imdata.deconvAlignedLick2 = struct;

% collect deconvolved data and lick data during protocol (collect trials deconvolved data)
LickStartSamplesPre = find(diff(sData.behavior.lickDs)==1)+1; % collect which sample licks start
LickStartSamplesPre(LickStartSamplesPre <= DataShownBeforeStimSamples) = NaN; % don't test the very beginning of recording
LickStartSamplesPre(LickStartSamplesPre > LickStartSamplesPre(end) - DataShownAfterStimSamples) = NaN;  % don't test the very end of recording
LickStartSamples = LickStartSamplesPre(~isnan(LickStartSamplesPre));
nLicks = length(LickStartSamples);

% run a moving mean on data to discard noise little bit 
Data = NaN(size(sData.imdata.roiSignals(2).dff));
for i=1:1:nROIs
    Data(i,:) = movmean(sData.imdata.roiSignals(2).dff(i,:),10);
end

% collect those signals when lick starts 
GrandCollection = NaN(nLicks-1*nROIs,DataShownLength);
for i=1:1:nROIs
    DataLicks{i} = NaN(nLicks-1,DataShownLength);
    for j=1:1:nLicks-1
        DataLicks{i}(j,1:DataShownLength) = Data(i,LickStartSamples(j)-DataShownBeforeStimSamples:LickStartSamples(j)+DataShownAfterStimSamples);
        GrandCollection((i*(nLicks-1)-(nLicks-1)+j),:) = Data(i,LickStartSamples(j)-DataShownBeforeStimSamples:LickStartSamples(j)+DataShownAfterStimSamples);
    end
end
DataLicks{nROIs + 1} = GrandCollection;

%%% PLOTS
% individual ROIS with all stimulation average data
Xstep = (DataShownAfterStimMs + DataShownBeforeStimMs) / (DataShownBeforeStimSamples+DataShownAfterStimSamples);
Xaxis = -Xstep*DataShownBeforeStimSamples:Xstep:Xstep*DataShownAfterStimSamples;

%{
%%% PLOT individual ROI activity
for roi=1:1:nROIs
    Ymax = max(max(DeconvDataLicks{roi}));
    figure('Color','white','visible',FigVisible)
    for i = 1:1:nLicks
        plot(Xaxis,DeconvDataLicks{roi}(i,:),'LineWidth',1.5); hold on;
    end
    line([0 0],[0 Ymax],'Color','black','LineStyle','--','LineWidth',1); hold on; 
    xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('activity (deconvolved signal)');
    %legend(legends{1},legends{2},legends{3},legends{4},'Location','northwest');  %legend{4},legend{5}, if more than 3 protocol exists
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi),' Lick-aligned activity'));
    FileName = strcat(sprintf('ROI-%d',roi),'-ROIAct-DeconvAlignedLick'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToLicks\ROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToLicks\ROIAct'),[FileName '.jpg'])));
    
end
close all;
%}

meanAllROIAlignedAct = NaN(nROIs,DataShownLength);
%%% PLOT overall mean ROI activity each porotocol in different plot  
for roi=1:1:nROIs
    meanAllROIAlignedAct(roi,:) = nanmean(DataLicks{roi},1);
    figure('Color','white','visible',FigVisible)
    plot(Xaxis,nanmean(DataLicks{roi},1),'LineWidth',1.5); hold on; % mean
    Ymax = max(nanmean(DataLicks{roi},1));
    line([0 0],[0 Ymax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
    xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('activity (dF/F)');
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi),' Mean lick-aligned activity'));
    FileName = strcat(sData.sessionInfo.fileID,'-Mean-lick-aligned-act',num2str(roi)); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),[FileName '.jpg'])));
end
close all
DataLicks{nROIs+2} = meanAllROIAlignedAct;

% Save file to same path where other files can be found 
save(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),strcat(sData.sessionInfo.fileID,'_dffToLicks.mat')),'DataLicks');

%%% plot mean of ROI's mean signal
figure('Color','white','visible',FigVisible)
plot(Xaxis,nanmean(meanAllROIAlignedAct,1),'Color','black','LineWidth',2); hold on; % mean
Ymax = max(max(meanAllROIAlignedAct));
line([0 0],[0 Ymax*1.1],'Color','red','LineStyle','--','LineWidth',1.5); hold on; 
hold on
for roi=1:1:nROIs
    plot(Xaxis,meanAllROIAlignedAct(roi,:),'Color',[0.75 0.75 0.75],'LineWidth',1); hold on; % mean
end
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity aligned to lick-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean lick-aligned activity all ROIs'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-lick-aligned-act-all-ROIs'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),[FileName '.jpg'])));


%%% plot all events' mean signal
figure('Color','white','visible',FigVisible)
plot(Xaxis,nanmean(GrandCollection,1),'Color','black','LineWidth',2); hold on; % mean
hold on
plot(Xaxis,nanmean(GrandCollection,1)+std(GrandCollection,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
plot(Xaxis,nanmean(GrandCollection,1)-std(GrandCollection,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
Ymax = max(nanmean(GrandCollection,1)+std(GrandCollection,1));
Ymin = min(nanmean(GrandCollection,1)-std(GrandCollection,1));
line([0 0],[Ymin Ymax],'Color','red','LineStyle','--','LineWidth',1.5); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity aligned to lick-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean lick-aligned activity all events'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-lick-aligned-act-all-events'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToLicks\ROIAct'),[FileName '.jpg'])));

 
end
