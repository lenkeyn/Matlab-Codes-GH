
function sData = NoraCVRconversionNoraML_CaData(sData,VelMin,IsDeconv)

%%% SET PARAMETERS 
VelMin = 0.1;
IsDenconv = 0;

nBins = sData.behavior.trialMatrices.meta.binNumber;
nTrials = sData.behavior.trialMatrices.meta.nAllTrials;
%FrameRate = sData.behavior.meta.imagingSamplingRate;

sData.imdataN = struct;
sData.imdataN.meta = struct;
sData.imdataN.binned = struct;
sData.imdataN.binned.VelMin = VelMin;
%sData.imdataN.binned.rewardPos = gol;

[nROIs, nSamples] = size(sData.imdata.roiSignals.dff);
sData.imdataN.roiSignals(2).dff(1:nROIs,1:nSamples) = sData.imdata.roiSignals.dff;
sData.imdataN.roiSignals(2).ch = 'green';
sData.imdataN.roiSignals(1).ch = 'red';
sData.imdataN.nROIs = nROIs;
sData.imdataN.nSamples = nSamples;

sData.behavior.lickDs(nSamples) = 0;


% saving folder
mkdir(sData.sessionInfo.savePath,'Imaging');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct');
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-dFF');

%{
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
    dff_lowpassLight(i,:) = filter(Hd,sData.imdataN.roiSignals(2).dff(i,1:nSamples)); 
end
sData.imdataN.roiSignals(2).dff_LPlight = single(dff_lowpassLight);
% plot slightly filtered transients
nROIOnFig = 10;
plotdata = dff_lowpassLight(:,1:nSamples);
% plot 10 ROIs to one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS
close all;
%}

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% binning
%sData.behavior.binning.samplesInBinIndex


% BINNING 
% calculating which sample belong to which bin
BinSize = sData.behavior.trialMatrices.meta.binSize;
TRNu = nTrials;
SampleInBin = NaN(nSamples,nBins);    
for j = 1:1:(nBins-1) % first search all data for the specified distance bin (indexes for entering a new bin can be calculated based on AbsDistDS data, j is distance cm, bin on the wheel (e.g. 0-1cm)
    for i = 1:1:nSamples % scan through the whole dataset and sort into bin j, between trials, there are NaNs on the column
       if  sData.behavior.signals.corridorPositionDs(i) >= (j-1)*BinSize && sData.behavior.signals.corridorPositionDs(i) < j*BinSize  %check ALL distanceDS datapoint if it belongs to the specified bin, if yes, continue
           SampleInBin(i,j) = i; % put sample index into the bin column. %MonIncr 
       end
    end 
end
if j == nBins-1 % last Bin, which is not a full bin or bigger, collect together what is before lapstart into last bin
   for i = 1:1:nSamples
      if  sData.behavior.signals.corridorPositionDs(i) > j*BinSize 
          SampleInBin(i,j+1) = i;          
      end
   end  
end
%%% correction for too fast trials, when one sample must be assign to two/more bins
% circularize data:
HeightPlus10Bin = size(SampleInBin,1)-1;
SampleInBinPlusBins = NaN(HeightPlus10Bin,nBins+10);
SampleInBinPlusBins(1:HeightPlus10Bin,1:nBins) = SampleInBin(1:HeightPlus10Bin,1:nBins);
SampleInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = SampleInBin(1:HeightPlus10Bin,1:10);
% correction
for j = 1:1:nSamples-2
    for i = 1:1:(nBins+10-1) 
       if  ~isnan(SampleInBinPlusBins(j,i)) && isnan(SampleInBinPlusBins(j+1,i))
           if isnan(SampleInBinPlusBins(j,i+1)) && isnan(SampleInBinPlusBins(j+1,i+1))
               SampleInBinPlusBins(j,i+1) = SampleInBinPlusBins(j,i); % if in a given bin no sample assign (too fast running), copy the previous bin's data there
           end
       end
    end
end
SampleInBinCorrect(1:HeightPlus10Bin,1:nBins) = SampleInBinPlusBins(1:HeightPlus10Bin,1:nBins);
SampleInBinCorrect(1:HeightPlus10Bin,1:5) = SampleInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+5);


% Make 3 matrices. EnterIntoBinSampleInd/LeaveBinSampleInd: when the animal enter/leave(in the next sample) a bin (sample ind, row: trial, col: bin); SampleSpentInBin: how many samples spend the aimal in this bin. At the end of recording there are NaNs in the matrices 
SampleInBinIsNaN1 = isnan(SampleInBinCorrect);
SampleInBinIsNaN2 = diff(SampleInBinIsNaN1,1)==-1; % matrix shows first sample in a bin. It is 1 and other is 0. (row: trial, column: bin)
if SampleInBinIsNaN1(1,1) == 0 % rare case when recording start in bin 1, it is needed to be changed in order to detect
    SampleInBinIsNaN2(1,1) = 1; 
end
SampleInBinIsNaN3 = diff(SampleInBinIsNaN1,1)==1; % matrix when last sample spent in a given bin, set to 1 and other is 0.
%LastSampleFix = find(SampleInBinIsNaN1(Samples,:) == 0);
LastSampleFix = SampleInBinIsNaN1(nSamples-1,:) == 0;  % I have to drop the last sample in order the code to function prperly
SampleInBinIsNaN3(nSamples-1,LastSampleFix) = 1 ; % drop last sample 
EnterIntoBinSampleInd = NaN(TRNu+1,nBins); % sample index when animal enter into given bin given trial
LeaveBinSampleInd = NaN(TRNu+1,nBins); % sample index when animal leaves a given bin
for i = 1:1:nBins % I need to do it in a complicated way using temporary arrays because if the last trial did not go to the end it gave error (mismatch in TRNu and bin-start at later bins)
    TempArray1 = zeros(TRNu+1,1);
    TempArray1 = find(SampleInBinIsNaN2(:,i)==1)+1; % first sample spend in a given bin, given trial
    if numel(TempArray1)<TRNu+1
        TempArray1(TRNu+1) = NaN;
    end
    EnterIntoBinSampleInd(:,i) = TempArray1;
    TempArray1 = zeros(TRNu+1,1);
    TempArray1 = find(SampleInBinIsNaN3(:,i)==1); % last sample spend in a given bin, given trial
    if numel(TempArray1)<TRNu+1
        TempArray1(TRNu+1) = NaN;
    end
    LeaveBinSampleInd(:,i) = TempArray1; % last sample spend in a given bin, given trial
end
SampleSpentInBin = LeaveBinSampleInd - EnterIntoBinSampleInd + 1; % spent in bin = last sample - first sample +1
% correction for bins in which the animal spent half bin time
% circularize matrix for calculation
HeightPlus10Bin = size(EnterIntoBinSampleInd,1)-1;
EnterIntoBinSampleIndPlusBins = NaN(HeightPlus10Bin,nBins+10);
EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,1:nBins) = EnterIntoBinSampleInd(1:HeightPlus10Bin,:);
EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = EnterIntoBinSampleInd(2:HeightPlus10Bin+1,1:10);
SampleSpentInBinPlusBins = NaN(HeightPlus10Bin,nBins+10);
SampleSpentInBinPlusBins(1:HeightPlus10Bin,1:nBins) = SampleSpentInBin(1:HeightPlus10Bin,:);
SampleSpentInBinPlusBins(1:HeightPlus10Bin,nBins+1:nBins+10) = SampleSpentInBin(2:HeightPlus10Bin+1,1:10);


for i = 1:1:TRNu
    for j = 1:1:nBins+10-1 %10 extra bin was added for circularization
        counter = 1;
        for k = 0:1:nBins+10-j-1
            if EnterIntoBinSampleIndPlusBins(i,j+k) == EnterIntoBinSampleIndPlusBins(i,j+k+1)
               counter = counter + 1;
            else
                break
            end
        end
        if counter > 1
            FrameSpent = SampleSpentInBinPlusBins(i,j)/counter;
            for k = 0:1:counter-1
                SampleSpentInBinPlusBins(i,j+k) = FrameSpent; % share the time spent there
            end
        end
    end
end
SampleSpentInBinCorrect(1:HeightPlus10Bin,1:nBins) = SampleSpentInBinPlusBins(1:HeightPlus10Bin,1:nBins);
SampleSpentInBinCorrect(2:HeightPlus10Bin,1:5) = SampleSpentInBinPlusBins(1:HeightPlus10Bin-1,nBins+1:nBins+5);

sData.behavior.binning.enterIntoBinIndex = EnterIntoBinSampleInd;
sData.behavior.binning.leaveBinIndex = LeaveBinSampleInd;
sData.behavior.binning.samplesSpentInBin = SampleSpentInBinCorrect;
sData.behavior.binning.enterIntoBinIndexExtended = EnterIntoBinSampleIndPlusBins;
sData.behavior.binning.SampleSpentInBinExtendedBins = SampleSpentInBinPlusBins; 
sData.behavior.binning.samplesInBinIndex = SampleInBinCorrect;




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
    MMSamplesLim{i} = NaN(MaxSamplesinBin,nBins); % for samples, creates i times a row-col matrices (row: max number of samples spent in a bin, column: Bin number)
    MMVelo{i} = NaN(MaxSamplesinBin,nBins); % same for all Velo data (wo limitation)
    MMVeloLim{i} = NaN(MaxSamplesinBin,nBins); % same for limited Velo data
end
% FILL UP MATRICES + Calculate average velocity within a bin and put into  CADATA.VeloInBin and CADATA.VeloLimInBin matrix
sData.imdataN.VeloInBin = NaN(nTrials-1,nBins); % mean real velo
sData.imdataN.VeloLimInBin = NaN(nTrials-1,nBins); % mean limited velo
for i = 1:1:nTrials-1
   for j = 1:1:nBins
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
      sData.imdataN.VeloLimInBin(i,j) = nanmean(MMVeloLim{i}(:,j)); % mean of limited velo
      sData.imdataN.VeloInBin(i,j) = nanmean(MMVelo{i}(:,j));  % mean of real velo
   end
end

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdataN.binned.RoidFF{i} = NaN(nTrials,nBins); 
    if IsDeconv >= 1
        sData.imdataN.binned.RoiDeconvolved{i} = NaN(nTrials,nBins);
        %sData.imdataN.binned.RoicumSpikes{i} = NaN(nTrials,nBin); 
        %sData.imdataN.binned.RoiiFreq{i} = NaN(nTrials,nBin); 
    end
end
 % load matrices with data in trial and bins
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:nBins
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
                   dFFInBin(m) = sData.imdataN.roiSignals(2).dff(i,(SampleInd+m-1)); % dFF
                   if IsDeconv >= 1
                      DeconvInBin(m) = sData.imdataN.roiSignals(2).deconv(i,(SampleInd+m-1)); % deconv 
                      %SpikeInBin(m) = sData.imdataN.roiSignals(2).cumSpikes(i,(SampleInd+m-1)); % spike rate
                      %iFreqInBin(m) = sData.imdataN.roiSignals(2).iFrequency(i,(SampleInd+m-1)); % 
                   end
                end
            end
            sData.imdataN.binned.RoidFF{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
            if IsDeconv >= 1
                sData.imdataN.binned.RoiDeconvolved{i}(j,k) = nanmean(DeconvInBin); % mean deconvolved Ca data 
                %sData.imdataN.binned.RoicumSpikes{i}(j,k) = nanmean(SpikeInBin); % mean spike rate of Ca data 
                %sData.imdataN.binned.RoiiFreq{i}(j,k) = nanmean(iFreqInBin);
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
sData.imdataN.roiStat = roiStat;

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% calculate and plot motion correction vectors for each bin
%sData = motionCorrVector(sData,filePath,savePath);

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
%nROIs = sData.imdataN.nROIs;
% plotdata: sData.imdataN.binned.RoidFF_SR_LP; sData.imdataN.binned.RoidFF_SR; sData.imdataN.binned.RoidFF

roiStart = 1;
roiEnd = nROIs;
FigVisible = 'off';

% plot binned dF/F 
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-dFF'); 
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdataN.binned.RoidFF{roi},1);

    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    plotHeatBinCa(sData.imdataN.binned.RoidFF{roi},sData.sessionInfo.fileID,roi,'dF/F',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdataN.binned.RoidFF{roi},1);
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
sData.imdataN.binned.MeanRoiAct = MeanRoiAct;
 
% mean signal RoiActBinned
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
figure('Color','white');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdataN.binned.MeanMeanRoiAct = nanmean(sData.imdataN.binned.MeanRoiAct);
Ymax = (max(sData.imdataN.binned.MeanMeanRoiAct(1,:)))*1.1;
Ymin = (min(sData.imdataN.binned.MeanMeanRoiAct(1,:)))*0.9;
plot(Xaxis,sData.imdataN.binned.MeanMeanRoiAct(1,1:nBins),'LineWidth',2)
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-dff');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

datatype = 0;
sData = placeCellMaosData2(sData,datatype);
% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end
