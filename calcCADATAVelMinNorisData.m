function sData = calcCADATAVelMinNorisData(sData,VelMin) 
%use this function to process Ca-imaging data 

%%% TO BE SET:  
% VelMin = 0.01; % minimum velocity. Below this, Ca-activity will be discarded

%%% recruiting Ca preprocessed data from ROIMANAGER: 
msgbox('Open dff file');
[~,filePath,~] = uigetfile('*.mat');
% dff = dff;

%%% SET PARAMETERS 
[nROIs, nSamples] = size(dff);
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
FrameRate = sData.behavior.meta.imagingSamplingRate;

sData.imdata = struct;
sData.imdata.binned = struct;
sData.imdata.meta = struct;

sData.imdata.nROIs = nROIs;
sData.imdata.nSamples = nSamples;

%%% load metadata
sData.imdata.meta = loadImagingMetadataNori;
load(fullfile(filePath,'_2p_stack_reg_ch1_001_rois.mat'));
sData.imdata.roiArray = struct;
sData.imdata.roiArray = roi_arr';  


%%% CONTROL. Compare frame signal in LV and number of frames in 2P recording:
if nSamples ~= sData.behavior.details.nFrames  
   msgbox(sprintf('Number of imaged frames and recorded by LV is not equal: %d vs %d',sData.behavior.details.nFrames,nSamples));
end
% set the smaller Sample-number for analysis
if nSamples >= sData.behavior.details.nFrames
   nSamples = sData.behavior.details.nFrames;
end

%%% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBin Matrix
SampleInBinVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end

% CREATE MATRIX for Velo-limited data based on velocity limit (discarded data will be represented by -1)
SampleInBinLim = sData.behavior.binning.samplesInBinIndex; % CADATA.SampleInBinLim = CADATA.SampleInBin;
% substitute indexes with low velo data with -1 : 
SampleInBinLim(sData.behavior.binning.veloBinned < VelMin) = -1; % CADATA.SampleInBinLim(CADATA.SampleInBinVelo < VelMin) = -1;

% PUT SAMPLE INDEXES (Velo limited) into multi-dimension matrix. 1st dim: trials, then columns are bins, rows contain samples in that trial and bin. 
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:)); % how many rows needed maximum in matrix (what was max time (samples) spend in a bin)   , MaxSamplesinBin = max(CADATA.SampleSpentInBin(:));
for i = 1:1:nTrials
    MMSamplesLim{i} = NaN(MaxSamplesinBin,BinNu); % for samples, creates i times a row-col matrices (row: max number of samples spent in a bin, column: Bin number)
    MMVelo{i} = NaN(MaxSamplesinBin,BinNu); % same for all Velo data (wo limitation)
    MMVeloLim{i} = NaN(MaxSamplesinBin,BinNu); % same for limited Velo data
end
% FILL UP MATRICES + Calculate average velocity within a bin and put into  CADATA.VeloInBin and CADATA.VeloLimInBin matrix
CADATA.VeloInBin = NaN (nTrials-1,BinNu); % mean real velo
CADATA.VeloLimInBin = NaN (nTrials-1,BinNu); % mean limited velo
for i = 1:1:nTrials-1
   for j = 1:1:BinNu
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
      CADATA.VeloLimInBin(i,j) = nanmean(MMVeloLim{i}(:,j)); % mean of limited velo
      CADATA.VeloInBin(i,j) = nanmean(MMVelo{i}(:,j));  % mean of real velo
   end
end


% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdata.binned.SingleRoidFF{i} = NaN(nTrials,BinNu); 
    sData.imdata.binned.SingleRoiDeconvolved{i} = NaN(nTrials,BinNu);
    sData.imdata.binned.SingleRoiSpikeRate{i} = NaN(nTrials,BinNu); 
end
 % load matrices with data in trial and bins
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:BinNu
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBin = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            dFFInBin = NaN(MaxSamplesinBin,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
            DeconvInBin = NaN(MaxSamplesinBin,1);
            SpikeInBin = NaN(MaxSamplesinBin,1);
            for m = 1:1:SpentInBin
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  dFFInBin(m) = CADATA.dFF(i,(SampleInd+m-1)); % dFF
                  DeconvInBin(m) = CADATA.Deconv(i,(SampleInd+m-1)); % deconv 
                  SpikeInBin(m) = CADATA.FiringRate(i,(SampleInd+m-1)); % spike rate
               end
            end
            sData.imdata.binned.SingleRoidFF{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
            sData.imdata.binned.SingleRoiDeconvolved{i}(j,k) = nanmean(DeconvInBin); % mean deconvolved Ca data 
            sData.imdata.binned.SingleRoiSpikeRate{i}(j,k) = nanmean(SpikeInBin); % mean spike rate of Ca data 
        end
    end
end

 
% PLOT HEATMAP of dff data in bins, plotHeatBinCa(data,fileID,roi,ylab,BinSize) 
% cue positions: %{
C1A = 23; C1B = 29; % velcro
C2A = 43; C2B = 49; % hot glue
C3A = 63; C3B = 69; % hot glue
C4A = 83; C4B = 89; % velcro

roiStart = 1;
roiEnd = nROIs;
% dF/F
for roi = roiStart:1:roiEnd
    plotHeatBinCa(sData.imdata.binned.SingleRoidFF{roi},sData.sessionInfo.fileID,roi,'dF/F',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--'); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--'); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--'); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--'); hold on;
    fname = strcat(sData.sessionInfo.fileID,'_roi',num2str(roi),'_dff');
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
end
close all;

%{
% deconvolved
for roi = roiStart:1:roiEnd
    plotHeatBin(sData.imdata.binned.SingleRoiDeconvolved{roi},CADATA.FilePath,roi,'Deconvolved signal',160,nTrials); 
    caxis([0 inf]); hold on;
    fname = strcat(num2str(roi),'deconv');
    savefig(fullfile(CADATA.FilePath,fname));
    saveas(gcf,(fullfile(CADATA.FilePath,[fname num2str(roi) '.jpg'])));
end
close all; %}


% spikes
for roi = roiStart:1:roiEnd
    plotHeatBinCa(sData.imdata.binned.SingleRoiSpikeRate{roi},sData.sessionInfo.fileID,roi,'Spike Rate',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on');
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--'); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--'); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--'); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--'); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--'); hold on;
    fname = strcat(sData.sessionInfo.fileID,'_roi',num2str(roi),'_spikerate');
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
end
close all;
%}

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end
