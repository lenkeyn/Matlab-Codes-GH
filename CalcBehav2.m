function sData = CalcBehav2(sData,nBins,IsOpto,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm) 

% calculate behavior data from TDMS data file and save it to behav file.  
% before this file use : sData.daqdata = loadTDMSdataNoriMateLV;

% set parameters:
% nBins = 30; % number of bins
% IsOptoSorted = 1; % was it an optical stimulated session? 0:no, 1:yes
% Sensitivity = 5; % sensitivity to detect optical stimulation
% If there was stimulation light in the (1/Sensitivity) percentage during the trial it should be considered as opto-on trial. 1/100 high sens, 1/5 low
% I had to introduce this parameter, since if I use full trial stimulation,the optical stimulation sometimes overlay into the beginning of next
% trial (for a few centimeters, until the animal receives the reward and licks). I do not want to recognize these trials as opto-on trials 
% (use sensitivity 5, and discard the beginning of the recordings). On the other hand when I am stimulating only a small part of the track (e.g.60-80 cm), 
% I want to detect these trials as opto-on, so I have to use % high sensitivity (100)
% SUM: set it = 5 for full trial stimulation, either square or sinus pulse (also if stimulation goes to the next lap, until reward), set to 100 to short stimulations - then it will be very sensitive


sData.behavior.meta.nBins = nBins;
sData.behavior.meta.binSize = 50*pi/nBins;        % double    Bin size when binning data. For instance, this can be the spatial bin size used for binning linear track data.

% Calculate Trial number (TRNu) and RewardStart array (array of zeros and ones, one when lap starts)
% Note: having 3 khz sampling, if the mouse run 100 cm/s, the photodiode signal (2cm black part) takes 20 ms, it is 60 samples. 
PhotoDiode = sData.daqdata.wheelDiode; % original analog data from tdms file
if any(PhotoDiode(:) > 2.0)
    PhotoDiode(PhotoDiode < 1.5) = 0; % binarized to zeros and ones
    PhotoDiode(PhotoDiode >= 1.5) = 1;
    PdZeroOne = diff(PhotoDiode) == 1; % search when PD signal change from zero to one - potential lap starts (and arteficial back-and-forth as well)
else
    PdZeroOne = PhotoDiode;
end
PotRewardStartIndex = find(PdZeroOne); % potential lap starts
RewardStartIndexPre = NaN(numel(PotRewardStartIndex),1); % real lap-starts
LapCounter = 0;
for i=2:1:numel(PotRewardStartIndex)  % the first 0-1 PD signal assumed to be real (but does not really matter if not)
    if sData.daqdata.distanceCm(PotRewardStartIndex(i)) - sData.daqdata.distanceCm(PotRewardStartIndex(i-1)) > 150 % lapstart if photodiode signal is 1 and animal runned at least 150 cm before previous reward
        LapCounter = LapCounter + 1; % sum up already runned laps
        RewardStartIndexPre(LapCounter) = PotRewardStartIndex(i);
    end
end
sData.behavior.details.rewardStartIndices = RewardStartIndexPre(~isnan(RewardStartIndexPre));
sData.behavior.wheelLap = LapCounter;

%frameMax = max(sData.daqdata.frameSignal);
%sData.daqdata.frameSignal(sData.daqdata.frameSignal > (frameMax/2)) = 1;
% Calculate number of imaging frames based on frame signal recorded by LabView
% generate fake frame signal (needed if there was no imaging but I want to downsample data the same way)
FrameStart = diff(sData.daqdata.frameSignal)== 1; % array with ones if the frame starts
if sum(FrameStart) == 0
    msgbox('No frame signal');
    FrameSignalDur = 98; %Samples
    nSamples = size(sData.daqdata.distanceCm,1);
    Samples = 1:1:nSamples;
    FakeFrameSignal = zeros(nSamples,1);
    FakeFrameSignal(rem(Samples,FrameSignalDur)==0) = 1;
    FrameStart = diff(FakeFrameSignal)== 1;
    FakeFrameSignal2 = 'yes';
end

sData.behavior.details.nFrames = sum(FrameStart); % frame number = sampleNumber usually behav.SampleNu = behav.FRNu;

% calculate number of frame signals and imaging frequency from frame signals  
FrameStartIndex = find(FrameStart); % array, indices when frame Starts , returns with the indices when FrameStart is not zero
sData.behavior.details.frameStartIndices = FrameStartIndex;
sData.behavior.meta.daqSamplingRate = sData.daqdata.meta.fs;
sData.behavior.meta.imagingSamplingRate = 1/((mean(diff(FrameStartIndex)))/sData.behavior.meta.daqSamplingRate); % frame frequency. Calculate mean sample number between frames, take the reciproc to get Hz

% cut recording before first reward and frame signals
Distance = sData.daqdata.distanceCm;
behavStart = min(sData.behavior.details.rewardStartIndices(sData.behavior.details.rewardStartIndices > FrameStartIndex(1)));
Distance(1:behavStart-1) = NaN;

% Aim: Calculate Absolute Distance (position on the wheel, zero at reward point) from cumulative distance data, downsaple to imaging sampling rate
% Calculate when rewards start, 
DistSubtract = Distance(sData.behavior.details.rewardStartIndices); % cumulative distances at lap starts (to be subtracted from cumDist data)
sData.behavior.stats.LapLengthCm = nanmean(diff(DistSubtract)); % mean length of laps, should be equal to Circumference
sData.behavior.stats.SDLapLengthCm = nanstd(diff(DistSubtract));
% inform the user
msgbox(sprintf('Full trial number is: %d. Mean LapLength: %g, Stdev LapLength: %g. Number of frames catched by LV frame signal: %g. Imaging sampling frequency is %g Hz.',sData.behavior.wheelLap,sData.behavior.stats.LapLengthCm,sData.behavior.stats.SDLapLengthCm,sData.behavior.details.nFrames,sData.behavior.meta.imagingSamplingRate));

% Calculate position on the wheel, zero in each lap is photodiode signal.
sData.behavior.wheelPos = NaN(numel(sData.daqdata.distanceCm),1);
LapCounter = 1;
for i = sData.behavior.details.rewardStartIndices(1):1:numel(sData.daqdata.distanceCm)
    if i > sData.behavior.details.rewardStartIndices(LapCounter+1)
        LapCounter = LapCounter + 1;
        if LapCounter == numel(sData.behavior.details.rewardStartIndices)
            break
        end
    end
    sData.behavior.wheelPos(i) = Distance(i) - DistSubtract(LapCounter);
end
% Downsample distance data to imaging frequency, consider only full laps with frame signal
sData.behavior.wheelPosDs = sData.behavior.wheelPos(FrameStartIndex);
% make it monotonic incr - in some cases,especially close to the reward, animals grab the wheel and move it back and forward, while their body is stationary. 
% Therefore in the position signal sometimes there are false backward movements. 
% I wanted to compensate for this and make the position signal monotonically increasing (if position moves backward, I gave the previous position value)
sData.behavior.wheelPosDsMonIncr = sData.behavior.wheelPosDs; % make it monotonically increasing 
for i=1:1:numel(sData.behavior.wheelPosDsMonIncr)-1
   if sData.behavior.wheelPosDsMonIncr(i)>150 && sData.behavior.wheelPosDsMonIncr(i+1)<5 % neglect lap start region
      continue
   elseif  sData.behavior.wheelPosDsMonIncr(i)> sData.behavior.wheelPosDsMonIncr(i+1) 
       sData.behavior.wheelPosDsMonIncr(i+1) = sData.behavior.wheelPosDsMonIncr(i);
   elseif  sData.behavior.wheelPosDsMonIncr(i) < 5   && sData.behavior.wheelPosDsMonIncr(i+1) >150 
       sData.behavior.wheelPosDsMonIncr(i+1) = sData.behavior.wheelPosDsMonIncr(i);
   end
end
% if the maximum backtracking is larger then 5 cm, discard the session
sData.behavior.wheelPosBackTrack = sData.behavior.wheelPosDsMonIncr - sData.behavior.wheelPosDs;
sData.behavior.details.MaxBackTRack = max(sData.behavior.wheelPosBackTrack);
msgbox(sprintf('The maximum backtracing was (cm): %g ',sData.behavior.details.MaxBackTRack));


% Velocity calculation
Circumfer = pi*50; % Wheel circumference
sData.behavior.stats.TheoreticLapLengthCm = Circumfer;
MovMean = 7 ; % moving window for smoothing
sData.behavior.runSpeedPre = movmean(diff(sData.daqdata.distanceCm),(sData.daqdata.meta.fs/5)); % smoothing for 3000/5=600 samples, 200 ms;
sData.behavior.runSpeed = sData.behavior.runSpeedPre / (1/sData.daqdata.meta.fs);


VelSubtr = diff(sData.behavior.wheelPosDsMonIncr); % calculate distance elapsed between frames
VelSubtr(1) = sData.behavior.wheelPosDsMonIncr(2); % first should have been calculated manually
% at lap starts have to adjust do not have negative valuse (157 -> 0)
for i=1:1:length(VelSubtr)
   if VelSubtr(i) < 0 && Circumfer > sData.behavior.wheelPosDsMonIncr(i)
       VelSubtr(i) = Circumfer - sData.behavior.wheelPosDsMonIncr(i) + sData.behavior.wheelPosDsMonIncr(i+1);
   elseif VelSubtr(i) < 0 % sometimes abs distance is calculated bigger than Circumference
       VelSubtr(i) = sData.behavior.wheelPosDsMonIncr(i+1);
   elseif VelSubtr(i) == 0
       VelSubtr(i) = 0.00001; % zeros is not good for later processing. 
   end
end
% the last value cannot be computed (now it is NAN), set it to the previous value, and also set a few more data to be able to calculate movmean
i = sum(~isnan(VelSubtr),1);  % search the last value which is non Nan
for j = 1:1:ceil(MovMean/2)
    VelSubtr(i+j) = VelSubtr(i);
end
VelPre = VelSubtr./(1/(sData.behavior.meta.imagingSamplingRate)); 
VelSmooth = movmean(VelPre,MovMean); % calculate movmean
NaNReplace = find(isnan(VelSmooth)); % in the beginning and end for movemean data, I have to substitute with original data
VelSmooth(NaNReplace,1) = VelPre(NaNReplace,1);
VelSmooth(VelSmooth<0) = 0.00001; % set negative values to almost zero
sData.behavior.runSpeedDs = VelSmooth;
sData.behavior.stats.VelMaxCmS = max(sData.behavior.runSpeedDs);

% generate downsampled lick data and downsampled water-given data
sData.behavior.lickDs = NaN(sData.behavior.details.nFrames-1,1); 
sData.behavior.waterRewardDs = NaN(sData.behavior.details.nFrames-1,1);
for i = 1:1:sData.behavior.details.nFrames-1
    TempLick = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for lick signal during a frame (scan)
    TempWater = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for water signal during a frame (scan)
    TempLick = sData.daqdata.lickSignal(FrameStartIndex(i):FrameStartIndex(i+1));
    TempWater = sData.daqdata.waterValve(FrameStartIndex(i):FrameStartIndex(i+1));
    if sum(TempLick)> max(diff(FrameStartIndex))/2 % if there is lick within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        sData.behavior.lickDs(i) = 1;
    else
        sData.behavior.lickDs(i) = 0;
    end
    if sum(TempWater) > max(diff(FrameStartIndex))/2 % if there is water within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        sData.behavior.waterRewardDs(i) = 1;
    else
        sData.behavior.waterRewardDs(i) = 0;
    end
end
sData.behavior.stats.LicksDSNu = sum(sData.behavior.lickDs);
sData.behavior.stats.WaterRewardNu = sum(sData.behavior.waterRewardDs);
% Check licking duarion
sData.behavior.stats.OneFrameDurInSample = 98;
LickStartTemp = find(diff(sData.daqdata.lickSignal)==1); % search lick signal start
LickEndTemp = find(diff(sData.daqdata.lickSignal)==-1); % lick signal ends
if LickEndTemp(1)> LickStartTemp(1)   % do not start during a lick
    LickEnd = LickEndTemp; % lick signal ends
else
    LickEnd = LickEndTemp(2:end);
end
if LickEnd(end)> LickStartTemp(end)   % do not end during a lick
    LickStart = LickStartTemp; % lick signal ends
else
    LickStart = LickStartTemp(1:end-1);
end
LickLengthSmpl = LickEnd - LickStart; % duration of lick signal
LickLengthSmpl = LickLengthSmpl(LickLengthSmpl> sData.behavior.stats.OneFrameDurInSample/5); % many single ones, the bottom of distribution is about 100
sData.behavior.stats.medianLickLengthSample = median(LickLengthSmpl);
sData.behavior.stats.SDLickLengthSample = std(LickLengthSmpl);

% Time axis
sData.behavior.timeInSec = ((1:length(sData.behavior.wheelPosDsMonIncr)-1)/sData.behavior.meta.imagingSamplingRate)';
sData.behavior.stats.RecordingDurationSec = max(sData.behavior.timeInSec);
sData.behavior.stats.RecordingDurationMin = max(sData.behavior.timeInSec)/60;

%%% BINNING
sData = binningNora(sData,nBins);
sData.behavior.wheelLapImaging = sum(sData.behavior.binning.enterIntoBinIndex(:,1)>0)-1;

%%% MAKING PLOTS
% save behavior plots into BEHAVIOR subfolder
mkdir(sData.sessionInfo.savePath,'Behavior');
savePath = strcat(sData.sessionInfo.savePath,'\Behavior');

% plot HeatBinVelo plot
[sData.behavior.binning.veloBinnedExtended, MeanVeloBinExtended,sData.behavior.WheelLapImagingExtended]= plotHeatBinVelo(sData,round(sData.behavior.stats.VelMaxCmS,-1)); % input= sData.behavior.runSpeedDs, Ymax for velo
FileName = strcat('VeloHeatBin-',sData.sessionInfo.fileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot HeatBinLicks plot
[sData.behavior.binning.lickBinnedExtended, MeanLickBinExtended,sData.behavior.binning.lickPerCmBinnedExtended,MeanLickCmExtended] = plotHeatBinLicksPlusBin(sData); 
FileName = strcat('LickHeatBin-',sData.sessionInfo.fileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot WaterGiven plot / hit-miss visible
SizeOfRewardZoneCm = 6;
[sData.behavior.binning.waterGivenBinnedExtended,sData.behavior.details.hitTrials,sData.behavior.details.hitRate] = plotHeatBinWaterGivenPlusBin(sData,SizeOfRewardZoneCm);
FileName = strcat('WaterGivenHeatBin-',sData.sessionInfo.fileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% write data into sData
if exist('FakeFrameSignal2')
    sData.behavior.meta.fakeframesignal = FakeFrameSignal2;
end
sData.behavior.binning.veloBinned = sData.behavior.binning.veloBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.lickBinned = sData.behavior.binning.lickBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.lickPerCmBinned = sData.behavior.binning.lickPerCmBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.waterGivenBinned = sData.behavior.binning.waterGivenBinnedExtended(1:sData.behavior.WheelLapImagingExtended,1:sData.behavior.meta.nBins);
sData.behavior.binning.meanVeloBinned = MeanVeloBinExtended(1:sData.behavior.meta.nBins);
sData.behavior.binning.meanLickBinned = MeanLickBinExtended(1:sData.behavior.meta.nBins);
sData.behavior.binning.meanLickPerCmBinned = MeanLickCmExtended(1:sData.behavior.meta.nBins);

%%% OPTICAL TRIAL SORTING 
% Plot stttings
RewardZoneLength = 6 ; % cm
Xstep = sData.behavior.meta.binSize;
XaxisEnd = nBins * sData.behavior.meta.binSize;
LickMax = 1.2*(max(MeanLickCmExtended(1:nBins)));
VMax = 1.2*(max(MeanVeloBinExtended(1:nBins)));


if IsOpto == 0 % if there was no optical stimulation in the session  
    %%% PLOT average lick / trials
    figure('Color','white'); 
    Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
    plot(Xaxis,MeanLickCmExtended(1:sData.behavior.meta.nBins)); hold on
    xlabel('Position on wheel (cm)');
    ylabel('Licks/cm');
    xlim([0 XaxisEnd])
    FileName = strcat('LickSum-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));  
    
    %%% PLOT average of trials speed
    figure('Color','white'); 
    Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
    plot(Xaxis,MeanVeloBinExtended(1:nBins)); hold on
    xlabel('Position on wheel (cm)');
    ylabel('Velocity (cm/s)');
    xlim([0 XaxisEnd])
    FileName = strcat('VeloSum-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
    
    if sum(sData.daqdata.optoSignal)>0
        msgbox(sprintf('Opto stimulation is detected, however session was set as non-optical session'));
    end
   
elseif IsOpto == 1 % if the session was optically stimulated

    % TRIAL SORTING 
    Opto = OptoTrialSorting3(sData,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm); 
    % OptoStimLimit = 9000; ( in msec): set the maximum optical stimulation duration (above this duration discard trials),DiscardCmBeginningCm = 10;  discard the first 10 cm for the opto sorting
    
    %%% PLOT average lick trial-types
    figure('Color','white'); 
    Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
    plot(Xaxis,Opto.MeanLick_OptoOff(1:nBins)); hold on;
    plot(Xaxis,Opto.MeanLick_OptoOn(1:nBins)); 
    plot(Xaxis,Opto.MeanLick_OptoAfter(1:nBins)); 
    xlabel('Position on wheel (cm)');
    ylabel('Licks/cm');
    xticks([0,25,50,75,100,125,150]);
    legend('Opto-Off','Opto-On','After-Opto','Location','north');
    FileName = strcat('LickSum-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    %%% PLOT Lick Heatmaps
    % OPTO-OFF
    plotdata = Opto.LickMatrix_OptoOff(1:(size(Opto.LickMatrix_OptoOff,1)),1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-Opto-Off'));
    FileName = strcat('LickHeat-OptoOff-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % OPTO-ON
    plotdata = Opto.LickMatrix_OptoOn(1:(size(Opto.LickMatrix_OptoOn,1)),1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-Opto-On'));
    FileName = strcat('LickHeat-OptoOn-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % AFTER-OPTO
    plotdata = Opto.LickMatrix_AfterOpto(1:(size(Opto.LickMatrix_AfterOpto,1)),1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-After-Opto'));
    FileName = strcat('LickHeat-AfterOpto-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


    %%%% SPEED

    %%% PLOT average of trials speed
    figure('Color','white'); 
    Xaxis = Xstep/2:Xstep:XaxisEnd-Xstep/2;
    plot(Xaxis,Opto.MeanVelo_OptoOff(1:nBins)); hold on;
    plot(Xaxis,Opto.MeanVelo_OptoOn(1:nBins)); hold on;
    plot(Xaxis,Opto.MeanVelo_OptoAfter(1:nBins)); hold on;
    xlabel('Position on wheel (cm)');
    ylabel('Velocity (cm/s)');
    xticks([0,25,50,75,100,125,150]);
    legend('Opto-Off','Opto-On','After-Opto','Location','South');
    FileName = strcat('VeloSum-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    %%% PLOT Velo Heatmaps
    % LIGHT-OFF
    plotdata = Opto.VeloMatrix_OptoOff(:,1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:nBins * sData.behavior.meta.binSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-Opto-Off'));
    FileName = strcat('VeloHeat-OptoOff-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % LIGHT-ON
    plotdata = Opto.VeloMatrix_OptoOn(:,1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:nBins * sData.behavior.meta.binSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-Opto-On'));
    FileName = strcat('VeloHeat-OptoOn-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % AFTER-OPTO
    plotdata = Opto.VeloMatrix_AfterOpto(:,1:nBins);
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:nBins * sData.behavior.meta.binSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID,'-After-Opto'));
    FileName = strcat('VeloHeat-AfterOpto-',sData.sessionInfo.fileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


    % write Optical stimulation data to sData
    sData.behavior.opto = Opto;
    sData.behavior.opto.IsOptoSession = 1;

end
  

% calculate behavioral performnace
if sData.behavior.wheelLap > 20
    sData = behavPerfomance(sData);
    sData = behavPerfomance2(sData);
end

% Save file to same path where LV files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end