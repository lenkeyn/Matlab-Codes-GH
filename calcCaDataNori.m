function sData = calcCaDataNori(sData,VelMin,laserPower,waveLength,fovCoord,pmtGainGreen,IsDeconv,gol) 
%use this function to process Ca-imaging data 

%%% TO BE SET:  
% VelMin = 0.01; % minimum velocity. Below this, Ca-activity will be discarded

%%% recruiting Ca preprocessed data from ROIMANAGER: 
msgbox('Open dff file');
[~,filePath,~] = uigetfile('*.mat');
% dff = dff;

%%% SET PARAMETERS 
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
FrameRate = sData.behavior.meta.imagingSamplingRate;

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

if IsDeconv == 1
    List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_cia_deconvolved.mat'));
    load(fullfile(filePath,List.name));
    sData.imdata.roiSignals(2).deconv = single(ciaDeconv); % cia_dec
end

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
List = dir(fullfile(filePath,'*_2p_stack_reg_ch1_001_signals.mat_cia_options.mat'));
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
%nROIOnFig = round((nROIs/3),-1)-1; % how many ROIs to see on the fig (best to set to X if rem(X/5)=1 , e.g. 19, 49, ...)
nROIOnFig = 10;
plotdata = dff_lowpassLight(:,1:nSamples);
% plot 10 ROIs to one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS

%{
% MEDIUM FILTERING to visualize very slow (minutes) transients
dff_lowpassMedium = NaN(nROIs,nSamples);
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.002,0.2,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
% Fp = 0.001 Shift = 500; Fp = 0.002 Shift = 200;   'Fp,Fst,Ap,Ast',0.002,0.2,1,60 SHift=200 ; 'Fp,Fst,Ap,Ast',0.0005,0.05,1,60 Shift=700
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    dff_lowpassMedium(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff(i,1:nSamples)); 
end
Shift = 200; % Shift the filtered data in X to aligned to real data
sData.imdata.roiSignals(2).dff_LPMedium = single(horzcat(dff_lowpassMedium(1:nROIs,Shift+1:nSamples),zeros(nROIs,Shift)));
plotdata = sData.imdata.roiSignals(2).dff_LPMedium;
%figure(); plot(dff_lowpassMedium(roi,1:nSamples-Shift)); hold on; plot(dff_lowpassMedium(roi,Shift:nSamples),'LineWidth',2);
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPMedium');


% STRONG FILTERING to visualize very slow (minutes) transients
dff_lowpassStrong = NaN(nROIs,nSamples);
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.0005,0.05,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
% Fp = 0.001 Shift = 500; Fp = 0.002 Shift = 200;   'Fp,Fst,Ap,Ast',0.002,0.2,1,60 SHift=200 ; 'Fp,Fst,Ap,Ast',0.0005,0.05,1,60 Shift=700
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    dff_lowpassStrong(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff(i,1:nSamples)); 
end
Shift = 700;
sData.imdata.roiSignals(2).dff_LPStrong = single(horzcat(dff_lowpassStrong(1:nROIs,Shift+1:nSamples),zeros(nROIs,Shift)));
plotdata = sData.imdata.roiSignals(2).dff_LPStrong;
%figure(); plot(dff_lowpassStrong(roi,1:nSamples-Shift)); hold on; plot(dff_lowpassStrong(roi,Shift:nSamples),'LineWidth',2);
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'ALL-NormROIsdFF-LPStrong');
%}

close all;


%%% Remove slow time-scale changes in fluoresence and baseline subtraction (5 percentile), baseline subtraction is not perfect (always above zero)
%Inspired Andreas by Dombeck et al (2010, Nat Neuro), code written by Andreas. Calculate basline in every 15 sec (in a moveing window) and subtract from data
RawNpilSubt = sData.imdata.roiSignals(2).roif - sData.imdata.roiSignals(2).npilf; 

Window1 = 5; % seconds, can be changed... , 10-15s, start: 10s
Percentile1 = 10; % I used 20%, but the transients became negative in baseline
Sampl = ceil(Window1*FrameRate); % samples to be used
SignalTemp1 = zeros(nROIs,nSamples); % temporary array for signal calculation
MeanBaseline = zeros(nROIs,1); % mean baseline for each ROI recording for calculate dFF
CollectBaseline = zeros(1,nSamples); % temporary collection of baseline value for a ROI
for i = 1:1:nROIs 
    % duplicate the beginning and end of the signal, and concatenate the first (and last) 15s of the signal to the beginning (and end) of the original signal
    SignalTemp2 = [RawNpilSubt(i,1:Sampl),RawNpilSubt(i,:),RawNpilSubt(i,(end-Sampl+1):end)]; % concatenate the duplicates and the original signal
    % Use a window of Window seconds around each data point to obtain X th percentile and subtract this from original signal to make the basline flat
    for j = Sampl:1:(nSamples+Sampl-1)
        Signal_window = SignalTemp2(1,(j-round(Sampl/2)):(j+round(Sampl/2))); % collect the data witihn the actual window to calculate baseline for datapoint 
        SignalTemp1(i,(j-Sampl+1)) = SignalTemp2(j) - prctile(Signal_window,Percentile1); % calculate X percentile of data and subtract
        CollectBaseline(1,(j-Sampl+1)) = prctile(Signal_window,Percentile1);
    end
    MeanBaseline(i,1) = mean(CollectBaseline(1,:)); % mean baseline for the whole ROI session, used for dF/F division
end
ROIsignals_raw_slow_removed = SignalTemp1; %

% Baseline already subtracted, generate dFF (devide with baseline)
dFF_slowRemoved = ROIsignals_raw_slow_removed ./ abs(MeanBaseline); % calculate dFF (baseline was subtracted in previous session)
sData.imdata.roiSignals(2).dff_slowRemoved = single(dFF_slowRemoved);

% calculate lowpassed data for visulization
dFF_slowRemoved_LP = NaN(nROIs,nSamples);
%Fp = 0.1; % PC neuron filtering for visualization
Fp = 0.02; % PV neurons strong filtering % test PV neurons (m8020-20190128-2): best to change only Fp between 0.1-0.01. Below 0.01 too much filtering, shift in time. Strong filering: Fp= 0.02, light filtering : Fp=0.1
d = fdesign.lowpass('Fp,Fst,Ap,Ast',Fp,1,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    dFF_slowRemoved_LP(i,:) = filter(Hd,dFF_slowRemoved(i,:)); 
end
sData.imdata.roiSignals(2).dff_slowRemoved_LP = single(dFF_slowRemoved_LP);

%nROIOnFig = round((nROIs/3),-1)-1; % how many ROIs to see on the fig (best to set to X if rem(X/5)=1 , e.g. 19, 49, ...)
plotdata = dFF_slowRemoved_LP;
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,savePath,'NormROIsdFF-slowRemLP');
close all;
%}

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

if IsDeconv == 2 
    %%% DECONVOLVE the calcium signals using fast non-negative deconvolution(?)
    % Be aware that deconvolving a signal that contains transients in signal due to movement or other artifacts may also make spikes in the deconvolved signal. 
    % Process deconvolveCa is developed by Pengcheng Zhou, Carnegie Mellon University, 2016. Ref: Friedrich J et.al., NIPS 2016, Fast Active Set Method for Online Spike Inference from Calcium Imaging. Also check Suite2p paper by Marius Pachitariu and Matteo Carandini, Kenneth D Harris https://doi.org/10.1101/061507
    %deconv_signal = NaN(nROIs,nSamples); % create emty temporary array
    for i = 1:1:nROIs
        SignalTemp3 = sData.imdata.roiSignals(2).dff_slowRemoved(i,1:nSamples); % temporary data container, use only selected ROIs
        SignalTemp3(SignalTemp3<0) = 0; % Set all negative values to zero
        [~,deconv] = deconvolveCa(double(SignalTemp3)); % use deconvolution script deconvolveCa. Without other input argument it will use 'constrained FOOPSI' method. It assumemes that the basline is zero(?). It calculated rise and decay time of transients and noise level as well for correct deconvolution, but you can set it manually.
        deconv_signal(i,:) = deconv;
    end
    sData.imdata.roiSignals(2).deconv = single(deconv_signal); %signal;
end

if IsDeconv >= 1
    %%% ESTIMATE 'SPIKING RATE', or cell activity (arbitrary 'firing rate') using the deconved signal
    % By integrating the deconvolved signal (it looks weird spiky reporting each activation) within a time window and setting a treshold what to consider activation, we can calculate some kind of arbitrary spike rate. At the end we apply gaussian smooting. Credit to Eivind Hennestad for most of this, and the remaining credit is for Andreas...
    % Set parmateres:
    CellFactor = 1; %multiplying factor for different cell tyes, PC = 1; VIP = 2;
    cumSamples = ceil(CellFactor*sData.imdata.meta.fps); % how many samples to use for 'memory' (cumulate signal) (31 frames = 1 sec). I set 0.5 sec (15 samples) for pyramidal cells, because rise of Ca-transient is usually within 0.5 sec. VIP: 2 sec
    %spikethresholdArray = max(sData.imdata.roiSignals(2).deconv,[],2); % treshold for consider an event as 'firing' (to discard small and fast fluorescent changes which does not seem real activation although the deconvolution method found them). Eivind uses 0.24, Andreas 0.11. I set 0.05, I found it reasonably good (maybe not detecting all single action potentials, but most of it does). Increase it to 0.24 might cause loosing 'few action potential' events, decreasing it to 0.01 might cause detection of a few noise deflections in the signal detected by deconvolution.
    spikethreshold = 0.1; %0.25
    % calculating spike rate
    estimated_firing_rate = zeros(nROIs,nSamples); % create empty array for data  
    for i = 1:1:nROIs
        %spikethreshold = spikethresholdArray(i)*0.2; 
        SignalTemp3 = sData.imdata.roiSignals(2).deconv(i,1:nSamples); % temporary array for calculating firing rate
        if isnan(SignalTemp3(i))
            continue
        end
        cumSum = 0; % The sum of integrating the spike rate
        lastSamples = zeros(1,round(cumSamples)); % Temporary array , the memory of last second. I want to reset the sum if nothing happened for the last second. Could maybe be shorter period.
        SignalTemp4 = zeros(size(SignalTemp3)); % Create vector to put spikes in
        for j = 1:1:length(SignalTemp3) % Run through the deconvolved signal
            lastSamples = horzcat(lastSamples(2:end),SignalTemp3(j)); % Put current sample at the end of "memory" vector
            if sum(lastSamples) == 0  % Integrate or reset.
                cumSum = 0;
            else 
                cumSum = cumSum + SignalTemp3(j);
            end
            if floor(cumSum/spikethreshold) >= 1 % Check if sum is over threshold. Count spikes and reset.
                SignalTemp4(j) = floor(cumSum/spikethreshold); %  SignalTemp4(j) = floor(cumSum/spikethreshold); 
                cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold; % cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold;
            end
            %{
            if floor(cumSum/spikethreshold) >= 1 % Check if sum is over threshold. Count spikes and reset.
                SignalTemp4(j) = floor(cumSum/spikethreshold); %  SignalTemp4(j) = floor(cumSum/spikethreshold); 
                cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold; % cumSum = cumSum - floor(cumSum/spikethreshold)*spikethreshold;
            end
            %}
        end
        estimated_firing_rate(i,:) = smoothdata(SignalTemp4,'gaussian',round(FrameRate/2));  %PC: round(FrameRate/2) % old:smoothdata(SignalTemp4./(1/FrameRate),'gaussian',round(FrameRate));smoothdata(SignalTemp4./(1/FrameRate),'gaussian',round(FrameRate));
    end
    sData.imdata.roiSignals(2).firingRate = single(estimated_firing_rate);
end

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
    sData.imdata.binned.RoidFF_SR{i} = NaN(nTrials,nBin);
    %sData.imdata.binned.RoidFF_SR_LP{i} = NaN(nTrials,nBin);
    if IsDeconv >= 1
    sData.imdata.binned.RoiDeconvolved{i} = NaN(nTrials,nBin);
    sData.imdata.binned.RoiSpikeRate{i} = NaN(nTrials,nBin); 
    %sData.imdata.binned.RoiDenoised{i} = NaN(nTrials,nBin); 
    %sData.imdata.binned.RoiDenoised_smoothed{i} = NaN(nTrials,nBin); 
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
            %dFFSRLPInBin = NaN(MaxSamplesinBin,1);
            if IsDeconv >= 1
            DeconvInBin = NaN(MaxSamplesinBin,1);
            %DenoisedInBin = NaN(MaxSamplesinBin,1);
            SpikeInBin = NaN(MaxSamplesinBin,1);
            %Denoised_smoothedInBin = NaN(MaxSamplesinBin,1);
            end
            for m = 1:1:SpentInBin
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  dFFInBin(m) = sData.imdata.roiSignals(2).dff(i,(SampleInd+m-1)); % dFF
                  dFFSRInBin(m) = sData.imdata.roiSignals(2).dff_slowRemoved(i,(SampleInd+m-1));
                  %dFFSRLPInBin(m) = sData.imdata.roiSignals(2).dff_slowRemoved_LP(i,(SampleInd+m-1)); % slow removed and filtered dFF
                  if IsDeconv >= 1
                  DeconvInBin(m) = sData.imdata.roiSignals(2).deconv(i,(SampleInd+m-1)); % deconv 
                  SpikeInBin(m) = sData.imdata.roiSignals(2).firingRate(i,(SampleInd+m-1)); % spike rate
                  %DenoisedInBin(m) = sData.imdata.roiSignals(2).denoised(i,(SampleInd+m-1)); % 
                  %Denoised_smoothedInBin(m) = sData.imdata.roiSignals(2).denoised_smoothed(i,(SampleInd+m-1)); % 
                  end
               end
            end
            sData.imdata.binned.RoidFF{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
            sData.imdata.binned.RoidFF_SR{i}(j,k) = nanmean(dFFSRInBin); %
            %sData.imdata.binned.RoidFF_SR_LP{i}(j,k) = nanmean(dFFSRLPInBin); %
            if IsDeconv >= 1
            sData.imdata.binned.RoiDeconvolved{i}(j,k) = nanmean(DeconvInBin); % mean deconvolved Ca data 
            sData.imdata.binned.RoiSpikeRate{i}(j,k) = nanmean(SpikeInBin); % mean spike rate of Ca data 
            %sData.imdata.binned.RoiDenoised{i}(j,k) = nanmean(DenoisedInBin); % mean spike rate of Ca data 
            %sData.imdata.binned.RoiDenoised_smoothed{i}(j,k) = nanmean(Denoised_smoothedInBin); % mean spike rate of Ca data 
            end
        end
    end
end

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
FigVisible = 'off';
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
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
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



% dF/F RoiActBinned_SR slow transients removed
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-SR'); % '\Imaging\RoiActBinned_SRLP'
MeanRoiActSR = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd
    MeanRoiActSR(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF_SR{roi},1);

    plotHeatBinCa(sData.imdata.binned.RoidFF_SR{roi},sData.sessionInfo.fileID,roi,'dF/F',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff-SR');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    MeanRoiActSR(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFF_SR{roi},1);
    Ymax = (max(MeanRoiActSR(roi,:)))*1.1;
    plot(Xaxis,MeanRoiActSR(roi,1:nBins),'LineWidth',2)
    line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-SR');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
   
end
close all;
sData.imdata.binned.MeanRoiAct_SR = MeanRoiActSR;
%}

%{
figure('Color','white')
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdata.binned.MeanMeanRoiAct_SR = nanmean(sData.imdata.binned.MeanRoiAct_SR);
Ymax = (max(sData.imdata.binned.MeanMeanRoiAct_SR(1,:)))*1.1;
Ymin = (min(sData.imdata.binned.MeanMeanRoiAct_SR(1,:)))*0.9;
plot(Xaxis,sData.imdata.binned.MeanMeanRoiAct_SR(1,1:nBins),'LineWidth',2)
line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs SR'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRoisSR-pos-tuning');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
%}


%%% plot smoothed firing rate (based on deconvolved signal)
if IsDeconv >= 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-deconv-FR');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-deconv-FR'); % '\Imaging\RoiActBinned_SRLP'
    MeanRoiActFR = NaN(nROIs,nBins);
    for roi = roiStart:1:roiEnd
        MeanRoiActFR(roi,1:nBins)= nanmean(sData.imdata.binned.RoiSpikeRate{roi},1);
    
        plotHeatBinCa(sData.imdata.binned.RoiSpikeRate{roi},sData.sessionInfo.fileID,roi,'firing rate',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
        caxis([0 inf]); hold on;
        line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-deconv');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        figure('Color','white','visible',FigVisible')
        Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
        %MeanRoiActFR(roi,1:nBins)= nanmean(sData.imdata.binned.RoiActFR{roi},1);
        Ymax = (max(MeanRoiActFR(roi,:)))+0.0001;
        plot(Xaxis,MeanRoiActFR(roi,1:nBins),'LineWidth',2)
        line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-deconv');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
        %}
    end
    close all;
    sData.imdata.binned.MeanRoiAct_FR = MeanRoiActFR;
    
    % mean pos tuning - FR
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    Ymax = (max(mean(MeanRoiActFR)))*1.1;
    Ymin = (min(mean(MeanRoiActFR)))*0.9;
    plot(Xaxis,mean(MeanRoiActFR),'LineWidth',2)
    line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs - FR'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-FR');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    
end


%}

%%% plot deconvolved signal
if IsDeconv >= 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-deconv');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-deconv'); % '\Imaging\RoiActBinned_SRLP'
    MeanRoiActDeconv = NaN(nROIs,nBins);
    for roi = roiStart:1:roiEnd
        MeanRoiActDeconv(roi,1:nBins)= nanmean(sData.imdata.binned.RoiDeconvolved{roi},1);
    
    %{
        if(any(isnan(MeanRoiActDeconv(roi,1:nBins))))
            continue 
        end
        plotHeatBinCa(sData.imdata.binned.RoiDeconvolved{roi},sData.sessionInfo.fileID,roi,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
        caxis([0 inf]); hold on;
        line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-deconv');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        figure('Color','white','visible',FigVisible')
        Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
        Ymax = (max(MeanRoiActDeconv(roi,:)))+0.0001;
        plot(Xaxis,MeanRoiActDeconv(roi,1:nBins),'LineWidth',2)
        line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-deconv');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
      %}
    end
    %}
    sData.imdata.binned.MeanRoiAct_Deconv = MeanRoiActDeconv;
    close all;
    
    % mean pos tuning
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    Ymax = (max(nanmean(MeanRoiActDeconv)))*1.1;
    Ymin = (min(nanmean(MeanRoiActDeconv)))*0.9;
    plot(Xaxis,nanmean(MeanRoiActDeconv),'LineWidth',2)
    line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs - deconvolved'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-deconv');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    
    % calculate all cells mean binned activity and plot in a heatplot 
    % convert cell array into 3d matrix, where x is row, y is column, z is ROI number
    sData.imdata.binned.RoiDeconvolvedMean = nanmean(cat(3,sData.imdata.binned.RoiDeconvolved{:}),3); 
    plotHeatBinCa(sData.imdata.binned.RoiDeconvolvedMean(1:nTrials,1:nBins),strcat(sData.sessionInfo.fileID,'-allROIsDeconv'),0,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    fname = strcat(sData.sessionInfo.fileID,'MeanDeconvBinned');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
end


%%% plot denoised signal
%{
if IsDeconv >= 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-denoised');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-denoised'); % '\Imaging\RoiActBinned_SRLP'
    MeanRoiActDenoised = NaN(nROIs,nBins);
    for roi = roiStart:1:roiEnd
        MeanRoiActDenoised(roi,1:nBins)= nanmean(sData.imdata.binned.RoiDenoised{roi},1);
    
        if(any(isnan(MeanRoiActDenoised(roi,1:nBins))))
            continue 
        end
        plotHeatBinCa(sData.imdata.binned.RoiDenoised{roi},sData.sessionInfo.fileID,roi,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
        caxis([0 inf]); hold on;
        line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-denoised');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        figure('Color','white','visible',FigVisible')
        Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
        Ymax = (max(MeanRoiActDenoised(roi,:)))+0.0001;
        plot(Xaxis,MeanRoiActDenoised(roi,1:nBins),'LineWidth',2)
        line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-denoised');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
      
    end
    sData.imdata.binned.MeanRoiAct_Denoised = MeanRoiActDenoised;
       
    close all;
    
    % mean pos tuning
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    Ymax = (max(nanmean(MeanRoiActDenoised)))*1.1;
    Ymin = (min(nanmean(MeanRoiActDenoised)))*0.9;
    plot(Xaxis,nanmean(MeanRoiActDenoised),'LineWidth',2)
    line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs - denoised'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-denoised');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    

    % calculate all cells mean binned activity and plot in a heatplot 
    % convert cell array into 3d matrix, where x is row, y is column, z is ROI number
    sData.imdata.binned.RoiDeconvolvedMean = nanmean(cat(3,sData.imdata.binned.RoiDeconvolved{:}),3); 
    plotHeatBinCa(sData.imdata.binned.RoiDeconvolvedMean(1:nTrials,1:nBins),strcat(sData.sessionInfo.fileID,'-allROIsDeconv'),0,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    fname = strcat(sData.sessionInfo.fileID,'MeanDeconvBinned');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
%}
%end

%%% plot denoised smoothed signal
%{
if IsDeconv >= 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-denoised-smoothed');
    savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-denoised-smoothed'); % '\Imaging\RoiActBinned_SRLP'
    MeanRoiActDenoisedSm = NaN(nROIs,nBins);
    for roi = roiStart:1:roiEnd
        MeanRoiActDenoisedSm(roi,1:nBins)= nanmean(sData.imdata.binned.RoiDenoised_smoothed{roi},1);
    
        if(any(isnan(MeanRoiActDenoisedSm(roi,1:nBins))))
            continue 
        end
        plotHeatBinCa(sData.imdata.binned.RoiDenoised_smoothed{roi},sData.sessionInfo.fileID,roi,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
        caxis([0 inf]); hold on;
        line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-denoisedSm');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        figure('Color','white','visible',FigVisible')
        Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
        Ymax = (max(MeanRoiActDenoisedSm(roi,:)))+0.0001;
        plot(Xaxis,MeanRoiActDenoisedSm(roi,1:nBins),'LineWidth',2)
        line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
        line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-denoisedSm');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
      
    end
    sData.imdata.binned.MeanRoiAct_Denoised_smoothed = MeanRoiActDenoisedSm;
    close all;
    
    
    % mean pos tuning
    figure();
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    Ymax = (max(nanmean(MeanRoiActDenoised)))*1.1;
    Ymin = (min(nanmean(MeanRoiActDenoised)))*0.9;
    plot(Xaxis,nanmean(MeanRoiActDenoised),'LineWidth',2)
    line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs - denoised'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-denoised');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    
    % calculate all cells mean binned activity and plot in a heatplot 
    % convert cell array into 3d matrix, where x is row, y is column, z is ROI number
    sData.imdata.binned.RoiDeconvolvedMean = nanmean(cat(3,sData.imdata.binned.RoiDeconvolved{:}),3); 
    plotHeatBinCa(sData.imdata.binned.RoiDeconvolvedMean(1:nTrials,1:nBins),strcat(sData.sessionInfo.fileID,'-allROIsDeconv'),0,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,'on'); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',1); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',1); hold on;
    fname = strcat(sData.sessionInfo.fileID,'MeanDeconvBinned');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
end
%}


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end

%{
A = max(sData.imdata.binned.MeanRoiAct,[],2);
B = sData.imdata.binned.MeanRoiAct./A;
plot(mean(B))
%}
