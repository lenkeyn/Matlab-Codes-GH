function sData  = CaDataAlignDeconvToStim(sData,FigVisible) 

%params:
DataShownAfterStimMs = 2000; % data shown after stimulus
DataShownBeforeStimMs = 2000; % data shown after stimulus 
DataShownAfterStimSamples = round(DataShownAfterStimMs/sData.imdata.meta.fps);
DataShownBeforeStimSamples = round(DataShownBeforeStimMs/sData.imdata.meta.fps);
DataShownLength = DataShownAfterStimSamples+DataShownBeforeStimSamples+1;

% set
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'AlignDeconvToStim');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),'ROIAct');
%SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim');
%nTrials = sData.behavior.wheelLapImaging-1;
nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
nOptoProt = max(sData.behavior.optoMoreProts.OptoStimProtTrials); %

% count the number of stimulated events in each protocol


%create struct
sData.imdata.deconvAlignedStim = struct;
realnOptoProt = NaN(nOptoProt,1); 

for pr = 2:1:nOptoProt %  repeat the process in all protocol types assuming that protocol 1 is the control
    sData.imdata.deconvAlignedStim.Protocol{pr} = struct;
    % collect deconvolved data and opto stim data during protocol (collect trials deconvolved data)
    ProtTrials = find(sData.behavior.optoMoreProts.OptoStimProtTrials == pr); % collect in which trials were the protocol used (first protocol is zero)
    ProtocolTrialStart = sData.behavior.binning.enterIntoBinIndex(ProtTrials,1);% collect in which samples a protocol (trial) starts (first column) and ends (second column), protocol changes in the beginning of trials
    ProtocolTrialEnd = sData.behavior.binning.leaveBinIndex(ProtTrials,nBins);
    % collect these signals and opto stimulus data 
    for i=1:1:nROIs
        counter = 1;
        for j=1:1:length(ProtTrials)
            DeconvDataStim(i,counter:counter+(ProtocolTrialEnd(j)-ProtocolTrialStart(j))) = sData.imdata.roiSignals(2).deconv(i,ProtocolTrialStart(j):ProtocolTrialEnd(j));
            if i == 1
                OptoStim(1,counter:counter+(ProtocolTrialEnd(j)-ProtocolTrialStart(j))) = sData.behavior.opto.LightOnSignalDS(ProtocolTrialStart(j):ProtocolTrialEnd(j));
                counter = counter + ProtocolTrialEnd(j)-ProtocolTrialStart(j);
            end
        end
    end
    % align deconv data to stimulus onset
    % look for stimulus onset:
    StimOnset = diff(OptoStim) == 1; %collects all stimulus start, but I want to find the first stimulus start in each trial:
    AllStims = sum(StimOnset);
    if AllStims == 0
        continue % if it was a control session, jump to next protocol
    end
    StimOnsetStart = NaN(length(ProtTrials),1);
    for i = 1:1:length(ProtTrials) % searh stimulus onset in each trial in this protocol
        TrialTemp =  sData.behavior.opto.LightOnSignalDS(ProtocolTrialStart(i):ProtocolTrialEnd(i));
        StimOnsetStart(i) = min(find(diff(TrialTemp)==1))+sData.behavior.binning.enterIntoBinIndex(ProtTrials(i),1);
    end
    nStims = length(ProtTrials); % if all individual stimulus should be taken into account then set it :  nStims = sum(StimOnset);
    realnOptoProt(pr,1) = pr; % how many protocol has opto stim
    %StimAlignData3D = NaN(nROIs,nStims,DataShownLength); % 3D matrix for storing aligned data
    %GrandAvrROIs = NaN(nROIs,DataShownLength);
    % load up data matrix
    for i = 1:1:nROIs
        %sData.imdata.deconvAlignedStim.ROIs{1,i} = NaN(nStims,DataShownLength);
        TempROIdata = NaN(nStims,DataShownLength);
        subtr = 0;
        for j = 1:1:nStims
            if StimOnsetStart(j)-DataShownBeforeStimSamples > 0  %length(DeconvDataStim) > (StimOnsetStart(j) + DataShownLength-1-subtr) && 
                %StimAlignData3D(i,j,1:DataShownLength) = sData.imdata.roiSignals(2).deconv(i,StimOnsetStart(j)-DataShownBeforeStimSamples:(StimOnsetStart(j) + DataShownBeforeStimSamples));
                TempROIdata(j,:) = sData.imdata.roiSignals(2).deconv(i,StimOnsetStart(j)-DataShownBeforeStimSamples+1:(StimOnsetStart(j)+DataShownBeforeStimSamples+1)); % individual ROI data for each stimulus
                TempROIdata(TempROIdata==Inf)=1;  % it is not a good solution, but I had to replace infinity numbers with something else....  
            end
            subtr = StimOnsetStart(j);
        end
        sData.imdata.deconvAlignedStim.Protocol{pr}.ROIData{i} = TempROIdata;
        sData.imdata.deconvAlignedStim.Protocol{pr}.GrandAvrROIs(i,:) = nanmean(TempROIdata,1);
        MaxValue{i}(1,pr) = max(nanmean(TempROIdata,1));
    end
    sData.imdata.deconvAlignedStim.Protocol{pr}.MeanOfAvrROIData = nanmean(sData.imdata.deconvAlignedStim.Protocol{pr}.GrandAvrROIs,1);
    sData.imdata.deconvAlignedStim.Protocol{pr}.StdOfAvrROIData = std(sData.imdata.deconvAlignedStim.Protocol{pr}.GrandAvrROIs,1);
    
    %clear DeconvDataStim OptoStim ProtTrials ProtocolTrialStart ProtocolTrialEnd StimOnset StimOnsetStart nStims
end

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

%%% PLOTS
% individual ROIS with all stimulation average data
Xstep = (DataShownAfterStimMs + DataShownBeforeStimMs) / (DataShownBeforeStimSamples+DataShownAfterStimSamples);
Xaxis = -Xstep*DataShownBeforeStimSamples:Xstep:Xstep*DataShownAfterStimSamples;
realOptoProt = realnOptoProt(~isnan(realnOptoProt));

%%% PLOT mean ROI activity in different protocols in one plot

for roi=1:1:nROIs
    Ymax = max(MaxValue{roi}(1,realOptoProt));
    figure('Color','white','visible',FigVisible)
    for i = 1:1:length(realOptoProt)
        plot(Xaxis,sData.imdata.deconvAlignedStim.Protocol{realOptoProt(i)}.GrandAvrROIs(roi,:),'LineWidth',1.5); hold on;
        legends{i} = strcat('Stim:',sData.stimProtocols(realOptoProt(i)).protocol);
    end
    line([0 0],[0 Ymax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
    xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('activity (deconvolved signal)');
    legend(legends{1},legends{2},legends{3},legends{4},'Location','northwest');  %legend{4},legend{5}, if more than 3 protocol exists
    title(strcat(sData.sessionInfo.fileID,'-ROI#',num2str(roi)));
    FileName = strcat(sprintf('ROI-%d',roi),'-ROIAct-DeconvAlignedStim'); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim\ROIAct'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim\ROIAct'),[FileName '.jpg'])));
    
end
close all;



%%% PLOT overall mean ROI activity each porotocol in different plot + sd 
for i = 1:1:length(realOptoProt)
    YMax = 0;
    YMin = max(max(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData));
    figure('Color','white','visible',FigVisible)
    plot(Xaxis,sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData,'LineWidth',1.5); hold on; % mean
    %{
    plot(Xaxis,(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData + sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData),'LineWidth',0.5,'Color','black'); % mean + sd
    plot(Xaxis,(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData - sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData),'LineWidth',0.5,'Color','black'); % mean - sd
    if YMax < max(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData + sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData)
        YMax = max(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData + sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData);
    end
    if YMin > min(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData - sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData)
       YMin = min(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData - sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.StdOfAvrROIData);
    end
    %}
    line([0 0],[YMin YMax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
    xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('activity (deconvolved signal)');
    legend(legends{i},'Location','northwest');  %legend{4},legend{5}, if more than 3 protocol exists
    title(strcat(sData.sessionInfo.fileID,'-Mean all ROIs data-Prot:',num2str(i)));
    FileName = strcat(sData.sessionInfo.fileID,'-AllROI-DeconvAlignedStim-Prot',num2str(i)); 
    savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),FileName));
    saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),[FileName '.jpg'])));
end

%%% PLOT overall mean ROI activity in different protocols in one plot (all ROIs in one plot)
YMax = 0;
figure('Color','white','visible',FigVisible)
for i = 1:1:length(realOptoProt)
    plot(Xaxis,sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData,'LineWidth',1.5); hold on;
    if YMax < max(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData)
        YMax = max(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData);
    end
end
line([0 0],[0 YMax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity (deconvolved signal)');
legend(legends{1},legends{2},legends{3},legends{4},'Location','northwest');  %legend{4},legend{5}, if more than 3 protocol exists
title(strcat(sData.sessionInfo.fileID,'-Mean all ROIs data'));
FileName = strcat(sData.sessionInfo.fileID,'-AllROI-DeconvAlignedStim'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),[FileName '.jpg'])));

%%% SAME, but smoothed data. PLOT overall mean ROI activity in different protocols in one plot (all ROIs in one plot)
YMax = 0;
Smoothing = 20;
figure('Color','white','visible',FigVisible)
for i = 1:1:length(realOptoProt)
    plot(Xaxis,smoothdata(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData,'gaussian',Smoothing),'LineWidth',1.5); hold on;
    if YMax < max(smoothdata(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData,'gaussian',Smoothing))
        YMax = max(smoothdata(sData.imdata.deconvAlignedStim.Protocol{1,realOptoProt(i)}.MeanOfAvrROIData,'gaussian',Smoothing));
    end
end
line([0 0],[0 YMax*1.1],'Color','black','LineStyle','--','LineWidth',1); hold on; 
xlabel('time (ms)'); ax = gca; ax.TickDir = 'out'; 
ylabel('activity (deconvolved signal)');
legend(legends{1},legends{2},legends{3},legends{4},'Location','northwest');  %legend{4},legend{5}, if more than 3 protocol exists
title(strcat(sData.sessionInfo.fileID,'-Mean all ROIs data-gauss-smoothed'));
FileName = strcat(sData.sessionInfo.fileID,'-AllROI-DeconvAlignedStim-GaussSmooth'); 
savefig(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),FileName));
saveas(gcf,(fullfile(strcat(sData.sessionInfo.savePath,'\Imaging\AlignDeconvToStim'),[FileName '.jpg'])));

 
end
