function sData = binningDeconvolved(sData) % signal: which signal to be binned: deconv, iFreq

% calculate binned data of 
signal = sData.imdata.roiSignals(2).deconv; 
str1 = 'deconv';

%mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),strcat('RoiActBinned-',str1));
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-',str1); 

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
nSamples = sData.imdata.meta.nFrames;
VelMin = 0.1;
 

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdata.binned.RoiDeconvolved{i} = NaN(nTrials,nBins);
end
MaxSamplesinBins = max(sData.behavior.binning.samplesSpentInBin(:));
%%% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBins Matrix
SampleInBinsVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinsVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end
SampleInBinsLim = sData.behavior.binning.samplesInBinIndex; 
SampleInBinsLim(SampleInBinsVelo < VelMin) = -1; 

 % load matrices with data in trial and bins
 % deconvolved and deconv +smoothed
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:nBins
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBins = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            deconvInBins = NaN(MaxSamplesinBins,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
            for m = 1:1:SpentInBins
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinsLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  deconvInBins(m) = signal(i,(SampleInd+m-1)); % dFF
                end
            end
            sData.imdata.binned.RoiDeconvolved{i}(j,k) = nanmean(deconvInBins); % mean of Ca data within a bin witihn a trial
        end
    end
end

% cue positions
C1A = sData.imdata.cues.C1A;
C1B = sData.imdata.cues.C1B;
C2A = sData.imdata.cues.C2A;
C2B = sData.imdata.cues.C2B;
C3A = sData.imdata.cues.C3A;
C3B = sData.imdata.cues.C3B;
C4A = sData.imdata.cues.C4A;
C4B = sData.imdata.cues.C4B;


roiStart = 1;
roiEnd = nROIs;
% plot data
FigVisible = 'off';
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoiDeconvolved{roi},1);
    %{
    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    plotHeatBinCa(sData.imdata.binned.RoiDeconvolved{roi},sData.sessionInfo.fileID,roi,'deconvolved activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-deconv');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoiDeconvolved{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1+0.0001;
    Ymax(Ymax==0) = 1;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    axis([0 160 0 Ymax]); % ceil(Ymax)
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-deconv');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
close all;
sData.imdata.binned.MeanRoiActDeconvolved = MeanRoiAct;

end