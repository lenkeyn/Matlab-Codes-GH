function sData  = CaDataAlignToReward(sData,FigVisible) 

FigVisible = 'on';
%params:
DataShownAfterStimMs = 1000; % data shown after stimulus
DataShownBeforeStimMs = 1000; % data shown before stimulus 
DataShownAfterStimSamples = round(DataShownAfterStimMs/sData.imdata.meta.fps);
DataShownBeforeStimSamples = round(DataShownBeforeStimMs/sData.imdata.meta.fps);
DataShownLength = DataShownAfterStimSamples+DataShownBeforeStimSamples+1;

% set
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'AlignDffToRew');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew'),'ROIAct');
%SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim');
%nTrials = sData.behavior.wheelLapImaging-1;
%nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;

%create struct
sData.imdata.AlignedRew = struct;

% collect deconvolved data and lick data during protocol (collect trials deconvolved data)
RewStartSamplesPre = find(diff(sData.behavior.waterRewardDs)==1)+1; % collect which sample licks start
RewStartSamplesPre(RewStartSamplesPre <= DataShownBeforeStimSamples) = NaN; % don't test the very beginning of recording
RewStartSamplesPre(RewStartSamplesPre > RewStartSamplesPre(end) - DataShownAfterStimSamples) = NaN;  % don't test the very end of recording
RewStartSamplesPre2 = RewStartSamplesPre;
for i = 2:1:length(RewStartSamplesPre)
    if RewStartSamplesPre(i) - RewStartSamplesPre(i-1) < DataShownBeforeStimSamples % there must be a reward-free period
        RewStartSamplesPre2(i) = NaN;
    end
end
RewStartSamples = RewStartSamplesPre2(~isnan(RewStartSamplesPre2));
nRew = length(RewStartSamples);

% run a moving mean on data to discard noise little bit 
Data = NaN(size(sData.imdata.roiSignals(2).dff));
for i=1:1:nROIs
    Data(i,:) = movmean(sData.imdata.roiSignals(2).dff(i,:),10);
end

% collect those signals when rew starts 
GrandCollection = NaN(nRew-1*nROIs,DataShownLength);
for i=1:1:nROIs
    DataRew{i} = NaN(nRew-1,DataShownLength);
    for j=1:1:nRew-1
        DataRew{i}(j,1:DataShownLength) = Data(i,RewStartSamples(j)-DataShownBeforeStimSamples:RewStartSamples(j)+DataShownAfterStimSamples);
        GrandCollection((i*(nRew-1)-(nRew-1)+j),:) = Data(i,RewStartSamples(j)-DataShownBeforeStimSamples:RewStartSamples(j)+DataShownAfterStimSamples);
    end
end
DataRew{nROIs + 1} = GrandCollection;

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
    meanAllROIAlignedAct(roi,:) = nanmean(DataRew{roi},1);
    figure('Color','white','visible',FigVisible)
    plot(Xaxis,nanmean(DataRew{roi},1),'LineWidth',1.5); hold on; % mean
    Ymax = max(nanmean(DataRew{roi},1));
    line([0 0],[0 Ymax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
    xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('activity (dF/F)');
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi),' Mean reward-aligned activity'));
    FileName = strcat(sData.sessionInfo.fileID,'-Mean-reward-aligned-act',num2str(roi)); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));
end
close all
DataRew{nROIs+2} = meanAllROIAlignedAct;

% Save file to same path where other files can be found 
save(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),strcat(sData.sessionInfo.fileID,'_dffToRew.mat')),'DataRew');

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
ylabel('activity aligned to reward-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean rew-aligned activity all ROIs'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-rew-aligned-act-all-ROIs'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));

%%% normalized ROI mean data
Max = max(meanAllROIAlignedAct,[],2);
meanAllROIAlignedActNorm = meanAllROIAlignedAct./Max;
% plot
figure('Color','white','visible',FigVisible)
for roi=1:1:nROIs
    plot(Xaxis,meanAllROIAlignedActNorm(roi,:),'Color',[0.75 0.75 0.75],'LineWidth',1); hold on; % mean
end
plot(Xaxis,nanmean(meanAllROIAlignedActNorm,1),'Color','black','LineWidth',2); hold on; % mean
Ymax = 1;
line([0 0],[0 Ymax],'Color','red','LineStyle','--','LineWidth',1.5); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity aligned to rew-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean rew-aligned activity all ROIs norm'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-rew-aligned-act-all-ROIs norm'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));



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
ylabel('activity aligned to reward-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean rew-aligned activity all events'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-rew-aligned-act-all-events'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));


% normalization of transients
Maxes = max(GrandCollection,[],2);
GrandCollectionNorm = GrandCollection./Maxes;
% plot
figure('Color','white','visible',FigVisible)
plot(Xaxis,nanmean(GrandCollectionNorm,1),'Color','black','LineWidth',2); hold on; % mean
hold on
plot(Xaxis,nanmean(GrandCollectionNorm,1)+std(GrandCollectionNorm,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
plot(Xaxis,nanmean(GrandCollectionNorm,1)-std(GrandCollectionNorm,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
Ymax = max(nanmean(GrandCollectionNorm,1)+std(GrandCollectionNorm,1));
Ymin = min(nanmean(GrandCollectionNorm,1)-std(GrandCollectionNorm,1));
line([0 0],[Ymin Ymax],'Color','red','LineStyle','--','LineWidth',1.5); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity aligned to reward-start (dF/F)');
title(strcat(sData.sessionInfo.fileID,' Mean reward-aligned normalized activity all events'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-rew-aligned-norm-act-all-events'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));


% speed aligned to reward
% run a moving mean on data to discard noise little bit 
SpeedData = NaN(1,sData.behavior.details.nFrames-1);
for i=1:1:nROIs
    SpeedData(i,:) = movmean(sData.behavior.runSpeedDs,10);
end
% claculate 
Speed = NaN(nRew-1,DataShownLength);
for j=1:1:nRew-1
    Speed(j,1:DataShownLength) = SpeedData(1,RewStartSamples(j)-DataShownBeforeStimSamples:RewStartSamples(j)+DataShownAfterStimSamples);
end

%%% PLOT mean speed   
figure('Color','white','visible',FigVisible)
plot(Xaxis,nanmean(Speed,1),'LineWidth',1.5); hold on; % mean
hold on
plot(Xaxis,nanmean(Speed,1)+nanstd(Speed,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
plot(Xaxis,nanmean(Speed,1)-nanstd(Speed,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
Ymax = max(nanmean(Speed,1)+nanstd(Speed,1));
line([0 0],[0 Ymax*1.1],'Color','red','LineStyle','--','LineWidth',1); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('speed (cm/s)');
title(strcat(sData.sessionInfo.fileID,' Average speed (cm/s)'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-reward-aligned-speed',num2str(roi)); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));


% licks aligned to reward
% run a moving mean on data to discard noise little bit 
LickData = NaN(1,sData.behavior.details.nFrames-1);
for i=1:1:nROIs
    LickData(i,:) = movmean(sData.behavior.lickDs,1);
end
% claculate 
Licks = NaN(nRew-1,DataShownLength);
for j=1:1:nRew-1
    Licks(j,1:DataShownLength) = LickData(1,RewStartSamples(j)-DataShownBeforeStimSamples:RewStartSamples(j)+DataShownAfterStimSamples);
end

%%% PLOT mean licks   
figure('Color','white','visible',FigVisible)
plot(Xaxis,nanmean(Licks,1),'LineWidth',1.5); hold on; % mean
hold on
plot(Xaxis,nanmean(Licks,1)+nanstd(Licks,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
plot(Xaxis,nanmean(Licks,1)-nanstd(Licks,1),'Color',[0.75 0.75 0.75],'LineStyle','- -','LineWidth',0.5); hold on; % mean
Ymax = max(nanmean(Licks,1)+nanstd(Licks,1));
line([0 0],[0 Ymax*1.1],'Color','red','LineStyle','--','LineWidth',1); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('licks');
title(strcat(sData.sessionInfo.fileID,' Average licks'));
FileName = strcat(sData.sessionInfo.fileID,'-Mean-reward-aligned-licks',num2str(roi)); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDffToRew\ROIAct'),[FileName '.jpg'])));
 
end
