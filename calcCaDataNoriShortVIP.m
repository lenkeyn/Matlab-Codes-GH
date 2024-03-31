function sData = calcCaDataNoriShortVIP(sData,VelMin,laserPower,waveLength,fovCoord,pmtGainGreen,IsDeconv,gol) 
%use this function to process Ca-imaging data 

%%% TO BE SET:  
% VelMin = 0.01; % minimum velocity. Below this, Ca-activity will be discarded
% IsDeconv = 1; % 2: do the deconvolution, 1: use ROI manager deconvolution, 0 : do not use deconvolution
% gol = 5; % gol#0 original reward, gol#2 50 cm forward, between hot glues, gol#3 50 cm before original after Velcro, gol = 10 no cues


%%% recruiting Ca preprocessed data from ROIMANAGER: 
msgbox('Open dff file');
[~,filePath,~] = uigetfile('*.mat');
% dff = dff;

%%% SET PARAMETERS 
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
%FrameRate = sData.behavior.meta.imagingSamplingRate;

sData.imdata = struct;
sData.imdata.meta = struct;
sData.imdata.binned = struct;
sData.imdata.binned.VelMin = VelMin;
sData.imdata.binned.rewardPos = gol;

%%% LOAD DATA
% load fluorescent raw data of ROIs
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_rois_meanf_raw.mat'));
load(fullfile(filePath,List.name)); %#ok<*LOAD>
sData.imdata.roiSignals(2).roif = single(roisMeanFRaw);
% load neuropil signal
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_npil_medif.mat'));
load(fullfile(filePath,List.name)); 
sData.imdata.roiSignals(2).npilf = single(npilMediF);
% load dff calculated by ROI manager
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_dff.mat'));
load(fullfile(filePath,List.name));
sData.imdata.roiSignals(2).dff = single(dff);
%sData.imdata.roiSignals(2).dff = single(ciaDenois);

%if IsDeconv == 1
    List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_cia_deconvolved.mat'));
    load(fullfile(filePath,List.name));
    sData.imdata.roiSignals(2).deconv = single(ciaDeconv); % cia_dec
%end

[nROIs, nSamples] = size(sData.imdata.roiSignals(2).dff);
sData.imdata.roiSignals(2).ch = 'green';
sData.imdata.roiSignals(1).ch = 'red';
sData.imdata.nROIs = nROIs;
sData.imdata.nSamples = nSamples;

% load signal extacion options
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_signal_extraction_options.mat'));
load(fullfile(filePath,List.name));
sData.imdata.signalExtractionOptions = options;
% load deconvolution options (if exists)
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_cia_options_ar1.mat'));
if ~isempty(List)
    load(fullfile(filePath,List.name));
    sData.imdata.signalExtractionOptions = ciaOptions{1,1};
end

% saving folder
mkdir(sData.sessionInfo.savePath,'Imaging');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-SR');

% load rois metadata
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_rois.mat'));
load(fullfile(filePath,List.name));
sData.imdata.roiArray = struct;
for i = 1:1:nROIs
    sData.imdata.roiArray(i).uid = roi_arr(1, i).uid;
    sData.imdata.roiArray(i).num = roi_arr(1, i).num;
    sData.imdata.roiArray(i).shape = roi_arr(1, i).shape;
    sData.imdata.roiArray(i).coordinates = roi_arr(1, i).coordinates;  
    sData.imdata.roiArray(i).imagesize = roi_arr(1, i).imagesize;
    sData.imdata.roiArray(i).center = roi_arr(1, i).center;
    sData.imdata.roiArray(i).area = roi_arr(1, i).area; 
    sData.imdata.roiArray(i).boundary = roi_arr(1, i).boundary;
    sData.imdata.roiArray(i).connectedrois = roi_arr(1, i).connectedrois;
    sData.imdata.roiArray(i).group = roi_arr(1, i).group;
    sData.imdata.roiArray(i).celltype = roi_arr(1, i).celltype;
    sData.imdata.roiArray(i).structure = roi_arr(1, i).structure;
    sData.imdata.roiArray(i).xyz = roi_arr(1, i).xyz;
    sData.imdata.roiArray(i).region = roi_arr(1, i).region;
    sData.imdata.roiArray(i).layer = roi_arr(1, i).layer;
    sData.imdata.roiArray(i).tags = roi_arr(1, i).tags;
    sData.imdata.roiArray(i).tag = roi_arr(1, i).tag;
    sData.imdata.roiArray(i).mask = roi_arr(1, i).mask;
end

% load metadata
sData.imdata.meta = loadImagingMetadataNori;
sData.imdata.meta.laserPower =  laserPower;          % mW                
sData.imdata.meta.waveLength =  waveLength;         % nm          
sData.imdata.meta.fovCoordinates = fovCoord;        % from center of window, AP ML (0 -500 means imaging at -2.2 mm left hemisphere)
sData.imdata.meta.pmtGain(1,1) = pmtGainGreen;

%%% CONTROL. Compare frame signal in LV and number of frames in 2P recording:
if nSamples ~= sData.behavior.details.nFrames  
   msgbox(sprintf('Number of recorded by LV frames and imaged frames is not equal: %d vs %d',sData.behavior.details.nFrames,nSamples));
end
% set the smaller Sample-number for analysis
if nSamples >= sData.behavior.details.nFrames
   nSamples = sData.behavior.details.nFrames;
end

% calculate slow removed transients, filtering Window in sec, percentile is for calculating baseline (ususlly 10-20%)
Window = 10; % based on m8056 VIP-Gcamp recording I find 10 sec good
Percentile = 20; % 20%
sData = SlowRemoved(sData,Window,Percentile);


%%% PLOT TRANSIENTS
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
% calculate lowpassed data for visulization
% LIGHT filtering to see better fast transients
dff_lowpassLight = NaN(nROIs,nSamples);
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,1,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    dff_lowpassLight(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff(i,1:nSamples)); 
    dff_lowpassLightSR(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff_slowRemoved(i,1:nSamples));
end
sData.imdata.roiSignals(2).dff_LPlight = single(dff_lowpassLight);
sData.imdata.roiSignals(2).dff_slowRemoved_LPlight = single(dff_lowpassLightSR);

% plot slightly filtered dff transients
%nROIOnFig = round((nROIs/3),-1)-1; % how many ROIs to see on the fig (best to set to X if rem(X/5)=1 , e.g. 19, 49, ...)
% plot 10 ROIs to one fig
nROIOnFig = 10; 
plotdata = dff_lowpassLight(:,1:nSamples);
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS
close all;


% plot slightly filtered dff_slowremoved transients
%nROIOnFig = round((nROIs/3),-1)-1; % how many ROIs to see on the fig (best to set to X if rem(X/5)=1 , e.g. 19, 49, ...)
% plot 10 ROIs to one fig
nROIOnFig = 10; 
plotdata = dff_lowpassLightSR(:,1:nSamples);
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'SR-NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'SR-ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS
close all;

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% calculating arbritary firing frequency
% calculates event weighted frequency of events that are separated by 0s
%{
if IsDeconv >= 1
    iFrequencyNorm = NaN(nROIs,nSamples); % create empty array for data  
    for i=287:1:nROIs
        signal = sData.imdata.roiSignals(2).deconv(i,1:nSamples);
        if max(signal) == 0
           continue 
        end
        iFrequency = nan(size(signal));
        eventPos = find(signal ~= 0);
        events = signal(eventPos);                          
        freques = sData.imdata.meta.fps./diff(eventPos);
        freques = [freques 0]; % add 0 as the last datapoint because frequency can only be calculated between two events 
        iFrequency(signal ~= 0) = freques.*events;          % fill data in the correct positions 
        iFrequency(1) = 0;                                  % add 0 as the first value for the fillmissing function
        iFrequency = fillmissing(iFrequency,'previous'); 
        iFrequencyNorm(i,1:nSamples) = iFrequency./min(signal(signal>0));
    end
    sData.imdata.roiSignals(2).iFrequencyNorm = single(iFrequencyNorm);
end
%}

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


%%% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBin Matrix
SampleInBinVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end

% CREATE MATRIX for Velo-limited data based on velocity limit (discarded data will be represented by -1)
SampleInBinLim = sData.behavior.binning.samplesInBinIndex; % CADATA.SampleInBinLim = CADATA.SampleInBin;
% substitute indexes with low velo data with -1 : 
SampleInBinLim(SampleInBinVelo < VelMin) = -1; % CADATA.SampleInBinLim(CADATA.SampleInBinVelo < VelMin) = -1;

% PUT SAMPLE INDEXES (Velo limited) into multi-dimension matrix. 1st dim: trials, then columns are bins, rows contain samples in that trial and bin. 
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:)); % how many rows needed maximum in matrix (what was max time (samples) spend in a bin)   , MaxSamplesinBin = max(CADATA.SampleSpentInBin(:));
for i = 1:1:nTrials
    MMSamplesLim{i} = NaN(MaxSamplesinBin,nBin); % for samples, creates i times a row-col matrices (row: max number of samples spent in a bin, column: Bin number)
    MMVelo{i} = NaN(MaxSamplesinBin,nBin); % same for all Velo data (wo limitation)
    MMVeloLim{i} = NaN(MaxSamplesinBin,nBin); % same for limited Velo data
end
% FILL UP MATRICES + Calculate average velocity within a bin and put into  CADATA.VeloInBin and CADATA.VeloLimInBin matrix
sData.imdata.VeloInBin = NaN(nTrials-1,nBin); % mean real velo
sData.imdata.VeloLimInBin = NaN(nTrials-1,nBin); % mean limited velo
for i = 1:1:nTrials-1
   for j = 1:1:nBin
      SampleInd = sData.behavior.binning.enterIntoBinIndex(i,j); % get sample index when enter into a bin  % SampleInd = CADATA.EnterIntoBinSampleInd(i,j); 
      if isnan(SampleInd)
         break 
      end
      SpentInBin = sData.behavior.binning.leaveBinIndex(i,j) - sData.behavior.binning.enterIntoBinIndex(i,j)+ 1 ; % get how many samples spent in that bin  % SpentInBin = CADATA.SampleSpentInBin(i,j);
      for k = 1:1:SpentInBin
          if SampleInBinLim(SampleInd+k-1,j) > 0 % limited values were set to -1, I do not want to contain them
             MMSamplesLim{i}(k,j) = SampleInd + k - 1; % limited samples
             MMVeloLim{i}(k,j) = SampleInBinVelo((SampleInd+k-1),j); % use only limited velo data  % CADATA.MMVeloLim{i}(k,j) = CADATA.SampleInBinVelo((SampleInd+k-1),j);
          end
          MMVelo{i}(k,j) = SampleInBinVelo((SampleInd+k-1),j); % use all velo data  % CADATA.MMVelo{i}(k,j) = CADATA.SampleInBinVelo((SampleInd+k-1),j); 
      end
      sData.imdata.VeloLimInBin(i,j) = nanmean(MMVeloLim{i}(:,j)); % mean of limited velo
      sData.imdata.VeloInBin(i,j) = nanmean(MMVelo{i}(:,j));  % mean of real velo
   end
end

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdata.binned.RoidFF{i} = NaN(nTrials,nBin);
    sData.imdata.binned.RoidFFSR{i} = NaN(nTrials,nBin);
    if IsDeconv >= 1
    sData.imdata.binned.RoiDeconvolved{i} = NaN(nTrials,nBin);
    sData.imdata.binned.RoiSpikeRate{i} = NaN(nTrials,nBin); 
    %sData.imdata.binned.iFrequencyNorm{i} = NaN(nTrials,nBin);
    end
end
 % load matrices with data in trial and bins
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:nBin
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBin = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            dFFInBin = NaN(MaxSamplesinBin,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
            dFFSRInBin = NaN(MaxSamplesinBin,1);
            if IsDeconv >= 1
            DeconvInBin = NaN(MaxSamplesinBin,1);
            %iFreqInBin = NaN(MaxSamplesinBin,1);
            end
            for m = 1:1:SpentInBin
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  dFFInBin(m) = sData.imdata.roiSignals(2).dff(i,(SampleInd+m-1)); % dFF
                  dFFSRInBin(m) = sData.imdata.roiSignals(2).dff_slowRemoved(i,(SampleInd+m-1));
                  if IsDeconv >= 1
                  DeconvInBin(m) = sData.imdata.roiSignals(2).deconv(i,(SampleInd+m-1)); % deconv 
                  %iFreqInBin(m) = sData.imdata.roiSignals(2).iFrequencyNorm(i,(SampleInd+m-1));
                  end
               end
            end
            sData.imdata.binned.RoidFF{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
            sData.imdata.binned.RoidFFSR{i}(j,k) = nanmean(dFFSRInBin); % mean of Ca data within a bin witihn a trial
            if IsDeconv >= 1
                sData.imdata.binned.RoiDeconvolved{i}(j,k) = nanmean(DeconvInBin); % mean deconvolved Ca data 
                %sData.imdata.binned.iFrequencyNorm{i}(j,k) = nanmean(iFreqInBin); 
            end
        end
    end
end
%}

% calculate ROIstats
roiStat = getRoiActivityStats(sData,2); % using channel 2
roiStat.meanPeakDff = nanmean(roiStat.peakDff);
roiStat.stdPeakDff = nanstd(roiStat.peakDff);
roiStat.meanSignalToNoise = nanmean(roiStat.signalToNoise);
roiStat.stdSsignalToNoise = nanstd(roiStat.signalToNoise);
roiStat.meanActivityLevel = nanmean(roiStat.activityLevel);
roiStat.stdActivityLevel = nanstd(roiStat.activityLevel);
sData.imdata.roiStat = roiStat;

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% calculate and plot motion correction vectors for each bin
sData = motionCorrVector(sData,filePath,savePath);

% PLOT HEATMAP of dff data in bins, plotHeatBinCa(data,fileID,roi,ylab,BinSize) 
% cue positions: %{
if gol == 0
C1A = 17; C1B = 23;  % velcro original set   C1A = 23; C1B = 29; Mate shifted reward in LV 6 cm: C1A = 17; C1B = 23; C1A = 10; C1B = 16;  %C1A = 15; C1B = 21; % MicheleLV era: C1A = 26; C1B = 36; 
C2A = 37; C2B = 43; % hot glue original set C2A = 43; C2B = 49; Mate shifted reward in LV 6 cm: C2A = 37; C2B = 43; C2A = 30; C2B = 36;  %C2A = 35; C2B = 41; % MicheleLV era: C2A = 58; C2B = 66;
C3A = 57; C3B = 63; % hot glue original set C3A = 63; C3B = 69; Mate shifted reward in LV 6 cm: C3A = 57; C3B = 63; C3A = 50; C3B = 56;  %C3A = 55; C3B = 61; % MicheleLV era: C3A = 86; C3B = 96;
C4A = 77; C4B = 83; % velcro original set   C4A = 83; C4B = 89; Mate shifted reward in LV 6 cm: C4A = 77; C4B = 83; C4A = 70; C4B = 76;  %C4A = 75; C4B = 81; % MicheleLV era: C4A = 119; C4B = 129;
elseif gol == 1
C1A = 67; C1B = 73;
C2A = 87; C2B = 93;
C3A = 107; C3B = 113;
C4A = 127; C4B = 133;
elseif gol == 2
C1A = 124; C1B = 130;
C2A = 144; C2B = 150;
C3A = 7; C3B = 13;
C4A = 27; C4B = 33;
%end
elseif gol == 4
%first cue setting: velcro, hot glue, hot glue, velcro
C1A = 26; C1B = 36;
C2A = 58; C2B = 66;
C3A = 88; C3B = 96;
C4A = 119; C4B = 129;
elseif gol == 5
% cue setting: sandpaper, hot glue, hot glue, sandpaper
C1A = 26; C1B = 26;
C2A = 46; C2B = 46;
C3A = 66; C3B = 66;
C4A = 86; C4B = 86;
elseif gol == 10
%no cues
C1A = 0; C1B = 0;
C2A = 0; C2B = 0;
C3A = 0; C3B = 0;
C4A = 0; C4B = 0;
end
%}
%C1A = 100; C1B = 102;

sData.imdata.cues.C1A = C1A;
sData.imdata.cues.C1B = C1B;
sData.imdata.cues.C2A = C2A;
sData.imdata.cues.C2B = C2B;
sData.imdata.cues.C3A = C3A;
sData.imdata.cues.C3B = C3B;
sData.imdata.cues.C4A = C4A;
sData.imdata.cues.C4B = C4B;

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
% plotdata: sData.imdata.binned.RoidFF_SR_LP; sData.imdata.binned.RoidFF_SR; sData.imdata.binned.RoidFF

roiStart = 1;
roiEnd = nROIs;
% dF/F RoiActBinned
FigVisible = 'on';
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned'); % '\Imaging\RoiActBinned_SRLP'
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF{roi},1);

    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    % sData.imdata.binned.RoidFF_SR_LP
    plotHeatBinCa(sData.imdata.binned.RoidFF{roi},sData.sessionInfo.fileID,roi,'dF/F',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    % sData.imdata.binned.RoidFF_SR_LP
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
close all;
sData.imdata.binned.MeanRoiAct = MeanRoiAct;


%%% PLot Mean position tuning of all ROIs
figure('Color','white');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdata.binned.MeanMeanRoiAct = nanmean(sData.imdata.binned.MeanRoiAct);
Ymax = (max(sData.imdata.binned.MeanMeanRoiAct(1,:)))*1.1;
Ymin = (min(sData.imdata.binned.MeanMeanRoiAct(1,:)))*0.9;
plot(Xaxis,sData.imdata.binned.MeanMeanRoiAct(1,1:nBins),'LineWidth',2)
line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

%%%Plot dFF_slowremoved signals
FigVisible = 'on';
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-SR'); % '\Imaging\RoiActBinned_SRLP'
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFFSR{roi},1);

    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    % sData.imdata.binned.RoidFF_SR
    plotHeatBinCa(sData.imdata.binned.RoidFFSR{roi},sData.sessionInfo.fileID,roi,'dF/F(SR)',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    % sData.imdata.binned.RoidFF_SR_LP
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFFSR{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
close all;
sData.imdata.binned.MeanRoiActSR = MeanRoiAct;


%%% PLot Mean position tuning of all ROIs slow removed transeints
figure('Color','white');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdata.binned.MeanMeanRoiActSR = nanmean(sData.imdata.binned.MeanRoiActSR);
Ymax = (max(sData.imdata.binned.MeanMeanRoiActSR(1,:)))*1.1;
Ymin = (min(sData.imdata.binned.MeanMeanRoiActSR(1,:)))*0.9;
plot(Xaxis,sData.imdata.binned.MeanMeanRoiActSR(1,1:nBins),'LineWidth',2)
line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs SR'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-SR');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));


%%%%
%{
%%% plot firing rate (based on deconvolved signal iFreq)
if IsDeconv >= 1
    FigVisible = 'on';
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-iFreq');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-iFreq'); 
    MeanRoiActiFreq = NaN(nROIs,nBins);
    for roi = roiStart:1:roiEnd
        MeanRoiActiFreq(roi,1:nBins)= nanmean(sData.imdata.binned.iFrequencyNorm{roi},1);
        plotHeatBinCa(sData.imdata.binned.iFrequencyNorm{roi},sData.sessionInfo.fileID,roi,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
        caxis([0 inf]); hold on;
        line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-iFreq');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        figure('Color','white','visible',FigVisible')
        Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
        Ymax = (max(MeanRoiActiFreq(roi,:)))+0.0001;
        plot(Xaxis,MeanRoiActiFreq(roi,1:nBins),'LineWidth',2)
        line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-iFreq');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
        
    end
    close all;
    sData.imdata.binned.MeanRoiAct_iFreq = MeanRoiActiFreq;
    
    % mean pos tuning - iFreq
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    Ymax = (max(mean(MeanRoiActiFreq)))*1.1;
    Ymin = (min(mean(MeanRoiActiFreq)))*0.9;
    plot(Xaxis,mean(MeanRoiActiFreq),'LineWidth',2)
    line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs - iFreq'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-iFreq');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    
end
%}
% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end


