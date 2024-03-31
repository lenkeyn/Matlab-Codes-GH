
function sData = NoraCconversionNoraML_CaData(sData,VelMin,IsDeconv)

%%% SET PARAMETERS 
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
%FrameRate = sData.behavior.meta.imagingSamplingRate;

sData.imdata = struct;
sData.imdata.meta = struct;
sData.imdata.binned = struct;
sData.imdata.binned.VelMin = VelMin;
%sData.imdata.binned.rewardPos = gol;

[nSamples, nROIs] = size(sData.RoiSignals.RoiSignals_Dff);
sData.imdata.roiSignals(2).dff(1:nROIs,1:nSamples) = sData.RoiSignals.RoiSignals_Dff';
sData.imdata.roiSignals(2).ch = 'green';
sData.imdata.roiSignals(1).ch = 'red';
sData.imdata.nROIs = nROIs;
sData.imdata.nSamples = nSamples;

sData.behavior.lickDs(nSamples) = 0;


% saving folder
mkdir(sData.sessionInfo.savePath,'Imaging');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-dFF');


%%% CONTROL. Compare frame signal in LV and number of frames in 2P recording:
if nSamples ~= sData.behavior.details.nFrames  
   msgbox(sprintf('Number of recorded by LV frames and imaged frames is not equal: %d vs %d',sData.behavior.details.nFrames,nSamples));
end
% set the smaller Sample-number for analysis
if nSamples >= sData.behavior.details.nFrames
   nSamples = sData.behavior.details.nFrames;
end


%%% PLOT TRANSIENTS
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');

% calculate lowpassed data for visulization
% LIGHT filtering to see better fast transients
dff_lowpassLight = NaN(nROIs,nSamples);
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,1,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    dff_lowpassLight(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff(i,1:nSamples)); 
end
sData.imdata.roiSignals(2).dff_LPlight = single(dff_lowpassLight);
% plot slightly filtered transients
nROIOnFig = 10;
plotdata = dff_lowpassLight(:,1:nSamples);
% plot 10 ROIs to one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS
close all;

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
    if IsDeconv >= 1
        sData.imdata.binned.RoiDeconvolved{i} = NaN(nTrials,nBin);
        %sData.imdata.binned.RoicumSpikes{i} = NaN(nTrials,nBin); 
        %sData.imdata.binned.RoiiFreq{i} = NaN(nTrials,nBin); 
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
            if IsDeconv >= 1
                DeconvInBin = NaN(MaxSamplesinBin,1);
                %iFreqInBin = NaN(MaxSamplesinBin,1);
                %SpikeInBin = NaN(MaxSamplesinBin,1);
            end
            for m = 1:1:SpentInBin
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                   dFFInBin(m) = sData.imdata.roiSignals(2).dff(i,(SampleInd+m-1)); % dFF
                   if IsDeconv >= 1
                      DeconvInBin(m) = sData.imdata.roiSignals(2).deconv(i,(SampleInd+m-1)); % deconv 
                      %SpikeInBin(m) = sData.imdata.roiSignals(2).cumSpikes(i,(SampleInd+m-1)); % spike rate
                      %iFreqInBin(m) = sData.imdata.roiSignals(2).iFrequency(i,(SampleInd+m-1)); % 
                   end
                end
            end
            sData.imdata.binned.RoidFF{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
            if IsDeconv >= 1
                sData.imdata.binned.RoiDeconvolved{i}(j,k) = nanmean(DeconvInBin); % mean deconvolved Ca data 
                %sData.imdata.binned.RoicumSpikes{i}(j,k) = nanmean(SpikeInBin); % mean spike rate of Ca data 
                %sData.imdata.binned.RoiiFreq{i}(j,k) = nanmean(iFreqInBin);
            end
        end
    end
end

% calculate ROIstats
roiStat = getRoiActivityStats(sData,2); % using channel 1 or 2
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
%sData = motionCorrVector(sData,filePath,savePath);

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
%nROIs = sData.imdata.nROIs;
% plotdata: sData.imdata.binned.RoidFF_SR_LP; sData.imdata.binned.RoidFF_SR; sData.imdata.binned.RoidFF

roiStart = 1;
roiEnd = nROIs;
FigVisible = 'off';

% plot binned dF/F 
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-dFF'); 
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF{roi},1);

    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    plotHeatBinCa(sData.imdata.binned.RoidFF{roi},sData.sessionInfo.fileID,roi,'dF/F',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1+0.0001;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
    if rem(roi,50)==0
        close all;
    end
end
close all;
sData.imdata.binned.MeanRoiAct = MeanRoiAct;
 
% mean signal RoiActBinned
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
figure('Color','white');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdata.binned.MeanMeanRoiAct = nanmean(sData.imdata.binned.MeanRoiAct);
Ymax = (max(sData.imdata.binned.MeanMeanRoiAct(1,:)))*1.1;
Ymin = (min(sData.imdata.binned.MeanMeanRoiAct(1,:)))*0.9;
plot(Xaxis,sData.imdata.binned.MeanMeanRoiAct(1,1:nBins),'LineWidth',2)
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-dff');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end
