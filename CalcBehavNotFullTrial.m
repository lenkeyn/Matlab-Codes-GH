function sData = CalcBehavNotFullTrial(sData,BinNu,IsOpto) 

% calculate LV data from TDMS data and save it to LVDATA file. From 2018.11.19. calculate lapstart when Water Valve opens (photodiode signal is not reliable)
% before this file use : sData.daqdata = loadTDMSdataNoriMateLV;
LapLengthSet = 90;

LVDATA = struct;
LVDATA.stats = struct;
LVDATA.BinNu = BinNu;
LVDATA.BinSize = LapLengthSet/BinNu;
%LVDATA.BinNu = ceil(50*pi/BinSize); % Bin Number,  circumference is 50*pi

% fileID
LVDATA.FileID = sData.sessionInfo.fileID;


% Calculate Trial number (TRNu) and RewardStart array (array of zeros and ones, one when lap starts)
% Note: having 3 khz sampling, if the mouse run 100 cm/s, the photodiode signal (2cm black part) takes 20 ms, it is 60 samples. 
PhotoDiode = sData.daqdata.wheelDiode; % original analog data from tdms file
if any(PhotoDiode(:) > 2.5)
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
LVDATA.PDStartIndex = RewardStartIndexPre(~isnan(RewardStartIndexPre));
%LVDATA.RewardStartIndex = RewardStartIndexPre(~isnan(RewardStartIndexPre));
%LVDATA.stats.TRNuFull = LapCounter;
%LVDATA.TRNu = LapCounter;
%}

% LapStart at WaterGiven
LVDATA.RewardStartIndex = find(diff(sData.daqdata.waterValve) == 1);
LVDATA.TRNu = sum(diff(sData.daqdata.waterValve) == 1)-2;
LVDATA.stats.TRNuFull = LVDATA.TRNu;


% Calculate number of frames = number of scans
%%% generate fake frame signal
FrameStart = diff(sData.daqdata.frameSignal)== 1; % array with ones if the frame starts
if sum(FrameStart) == 0
    msgbox('No frame signal');
    FrameSignalDur = 98; %Samples
    nSamples = size(sData.daqdata.distanceCm,1);
    Samples = 1:1:nSamples;
    FakeFrameSignal = zeros(nSamples,1);
    FakeFrameSignal(rem(Samples,FrameSignalDur)==0) = 1;
    FrameStart = diff(FakeFrameSignal)== 1;
    LVDATA.FakeFrameSignal = 'yes';
end

LVDATA.FrameNu = sum(FrameStart); % frame number = sampleNumber usually
%LVDATA.SampleNu = LVDATA.FRNu;

% calculate number of frame signals and imaging frequency from frame signals  
FrameStartIndex = find(FrameStart); % array, indices when frame Starts , returns with the indices when FrameStart is not zero
LVDATA.FrameStartIndex = FrameStartIndex;
LVDATA.DaqSamplingRate = sData.daqdata.meta.fs;
LVDATA.FrameSignalRate = 1/((mean(diff(FrameStartIndex)))/LVDATA.DaqSamplingRate); % frame frequency. Calculate mean sample number between frames, take the reciproc to get Hz

% cut LV recording before first reward and frame signals
Distance = sData.daqdata.distanceCm;
%Distance(1:FrameStartIndex(1)-1) = NaN;
Start = min(LVDATA.RewardStartIndex(LVDATA.RewardStartIndex > FrameStartIndex(1)));
Distance(1:Start-1) = NaN;

% Aim: Calculate Absolute Distance (zero at reward point) from cumulative distance data, downsaple to imaging sampling rate, FrameSignalRate
% Calculate when rewards start, 
DistSubtract = Distance(LVDATA.RewardStartIndex); % cumulative distances at lap starts (to be subtracted from cumDist data)
LVDATA.stats.LapLengthCm = nanmean(diff(DistSubtract)); % mean length of laps, should be equal to Circumference
LVDATA.stats.SDLapLengthCm = nanstd(diff(DistSubtract));
% inform the user
msgbox(sprintf('Full trial number is: %d. Mean LapLength: %g, Stdev LapLength: %g. Number of frames catched by LV frame signal: %g. Imaging sampling frequency is %g Hz.',LVDATA.stats.TRNuFull,LVDATA.stats.LapLengthCm,LVDATA.stats.SDLapLengthCm,LVDATA.FrameNu,LVDATA.FrameSignalRate));

% Calculate position on the wheel, zero in each lap is photodiode/watervalve opening signal start.
AbsDist = NaN(numel(sData.daqdata.distanceCm),1);
LapCounter = 1;
for i = LVDATA.RewardStartIndex(1):1:numel(sData.daqdata.distanceCm)
    if i > LVDATA.RewardStartIndex(LapCounter+1)
        LapCounter = LapCounter + 1;
        if LapCounter == numel(LVDATA.RewardStartIndex)
            break
        end
    end
    AbsDist(i) = Distance(i) - DistSubtract(LapCounter);
end
% Downsample distance dadta to imaging frequency, consider only full laps with frame signal
LVDATA.AbsDistDS = AbsDist(FrameStartIndex);
% make it monotonic incr
LVDATA.AbsDistDSMonIncr = LVDATA.AbsDistDS; % make it monotonically increasing 
for i=1:1:numel(LVDATA.AbsDistDSMonIncr)-1
   if LVDATA.AbsDistDSMonIncr(i)>LapLengthSet-5 && LVDATA.AbsDistDSMonIncr(i+1)<5
      continue
   elseif  LVDATA.AbsDistDSMonIncr(i)> LVDATA.AbsDistDSMonIncr(i+1) 
       LVDATA.AbsDistDSMonIncr(i+1) = LVDATA.AbsDistDSMonIncr(i);
   elseif  LVDATA.AbsDistDSMonIncr(i) < 5   && LVDATA.AbsDistDSMonIncr(i+1) >LapLengthSet-5 
       LVDATA.AbsDistDSMonIncr(i+1) = LVDATA.AbsDistDSMonIncr(i);
   end
end


% Velocity calculation
Circumfer = pi*50; % Wheel circumference
LVDATA.stats.TheoreticLapLengthCm = Circumfer;
MovMean = 7 ; % moving window for smoothing
%%% calculation nonDS Velocity data
VelocityInstPre = movmean(diff(sData.daqdata.distanceCm),(sData.daqdata.meta.fs/5)); % smoothing for 3000/5=600 samples, 200 ms;
VelocityInst = VelocityInstPre / (1/sData.daqdata.meta.fs);


VelSubtr = diff(LVDATA.AbsDistDSMonIncr); % calculate distance elapsed between frames
VelSubtr(1) = LVDATA.AbsDistDSMonIncr(2); % first should have been calculated manually
% at lap starts have to adjust do not have negative valuse (157 -> 0)
for i=1:1:length(VelSubtr)
   if VelSubtr(i) < 0 && Circumfer > LVDATA.AbsDistDSMonIncr(i)
       VelSubtr(i) = Circumfer - LVDATA.AbsDistDSMonIncr(i) + LVDATA.AbsDistDSMonIncr(i+1);
   elseif VelSubtr(i) < 0 % sometimes abs distance is calculated bigger than Circumference
       VelSubtr(i) = LVDATA.AbsDistDSMonIncr(i+1);
   elseif VelSubtr(i) == 0
       VelSubtr(i) = 0.00001; % zeros is not good for later processing. 
   end
end
% the last value cannot be computed (now it is NAN), set it to the previous value, and also set a few more data to be able to calculate movmean
i = sum(~isnan(VelSubtr),1);  % search the last value which is non Nan
for j = 1:1:ceil(MovMean/2)
    VelSubtr(i+j) = VelSubtr(i);
end
VelPre = VelSubtr./(1/(LVDATA.FrameSignalRate)); 
VelSmooth = movmean(VelPre,MovMean); % calculate movmean
NaNReplace = find(isnan(VelSmooth)); % in the beginning and end movemean loose data, I have to substitute with original data
VelSmooth(NaNReplace,1) = VelPre(NaNReplace,1);
VelSmooth(VelSmooth<0) = 0.00001; % set negative values to zero
LVDATA.VelDS = VelSmooth;
LVDATA.stats.VelMaxCmS = max(LVDATA.VelDS);

% generate downsampled lick data and downsampled water-given data
LVDATA.LickDS = NaN(LVDATA.FrameNu-1,1); 
LVDATA.WaterRewardDS = NaN(LVDATA.FrameNu-1,1);
for i = 1:1:LVDATA.FrameNu-1
    TempLick = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for lick signal during a frame (scan)
    TempWater = NaN(max(diff(FrameStartIndex))+1,1); % temporary array for water signal during a frame (scan)
    TempLick = sData.daqdata.lickSignal(FrameStartIndex(i):FrameStartIndex(i+1));
    TempWater = sData.daqdata.waterValve(FrameStartIndex(i):FrameStartIndex(i+1));
    if sum(TempLick)> mean(diff(FrameStartIndex))/2 % if there is lick within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        LVDATA.LickDS(i) = 1;
    else
        LVDATA.LickDS(i) = 0;
    end
    if sum(TempWater) > mean(diff(FrameStartIndex))/2 % if there is water within this frame, put 1 into DS lick array (if lick signal is one more than half time during frame scanning (one frame is 320 samples, one lick is usually 500 samples, sampling is 10 kHz)
        LVDATA.WaterRewardDS(i) = 1;
    else
        LVDATA.WaterRewardDS(i) = 0;
    end
end
LVDATA.stats.LicksDSNu = sum(LVDATA.LickDS);
LVDATA.stats.WaterRewardNu = sum(LVDATA.WaterRewardDS);
% Check licklength
LVDATA.stats.OneFrameDurInSample = 98;
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
LickLengthSmpl = LickEnd - LickStart; % length of lick signal
LickLengthSmpl = LickLengthSmpl(LickLengthSmpl> LVDATA.stats.OneFrameDurInSample/5); % many single ones, the bottom of distribution is about 100
LVDATA.stats.medianLickLengthSample = median(LickLengthSmpl);
LVDATA.stats.SDLickLengthSample = std(LickLengthSmpl);

% Time axis
LVDATA.TimeSec = ((1:length(LVDATA.AbsDistDSMonIncr)-1)/LVDATA.FrameSignalRate)';
LVDATA.stats.RecordingDurationSec = max(LVDATA.TimeSec);
LVDATA.stats.RecordingDurationMin = max(LVDATA.TimeSec)/60;

% BINNING
% calculating which sample belong to which bin
SampleInBin = NaN(LVDATA.FrameNu,LVDATA.BinNu);    
for j = 1:1:(LVDATA.BinNu-1) % first search all data for the specified distance bin (indexes for entering a new bin can be calculated based on AbsDistDS data, j is distance cm, bin on the wheel (e.g. 0-1cm)
    for i = 1:1:LVDATA.FrameNu % scan through the whole dataset and sort into bin j, between trials, there are NaNs on the column
       if  LVDATA.AbsDistDSMonIncr(i) >= (j-1)*LVDATA.BinSize && LVDATA.AbsDistDSMonIncr(i) < j*LVDATA.BinSize  %check ALL distanceDS datapoint if it belongs to the specified bin, if yes, continue
           SampleInBin(i,j) = i; % put sample index into the bin column. 
       end
    end 
end
if j == LVDATA.BinNu-1 % last Bin, which is not a full bin or bigger, collect together what is before lapstart into last bin
   for i = 1:1:LVDATA.FrameNu
      if  LVDATA.AbsDistDSMonIncr(i) > j*LVDATA.BinSize 
          SampleInBin(i,j+1) = i;          
      end
   end  
end
LVDATA.SampleInBinMatrix = SampleInBin;

% Make 3 matrices. EnterIntoBinSampleInd/LeaveBinSampleInd: when the animal enter/leave(in the next sample) a bin (sample ind, row: trial, col: bin); SampleSpentInBin: how many samples spend the aimal in this bin. At the end of recording there are NaNs in the matrices 
SampleInBinIsNaN1 = isnan(SampleInBin);
SampleInBinIsNaN2 = diff(SampleInBinIsNaN1,1)==-1; % matrix shows first sample in a bin. It is 1 and other is 0. (row: trial, column: bin)
if SampleInBinIsNaN1(1,1) == 0 % rare case when recording start in bin 1, it is needed to be changed in order to detect
    SampleInBinIsNaN2(1,1) = 1; 
end
SampleInBinIsNaN3 = diff(SampleInBinIsNaN1,1)==1; % matrix when last sample spent in a given bin, set to 1 and other is 0.
%LastSampleFix = find(SampleInBinIsNaN1(LVDATA.FrameNu,:) == 0);
LastSampleFix = SampleInBinIsNaN1(LVDATA.FrameNu,:) == 0;  % I have to drop the last sample in order the code to function prperly
SampleInBinIsNaN3(LVDATA.FrameNu-1,LastSampleFix) = 1 ; % drop last sample 
EnterIntoBinSampleInd = NaN(LVDATA.TRNu+1,LVDATA.BinNu); % sample index when animal enter into given bin given trial
LeaveBinSampleInd = NaN(LVDATA.TRNu+1,LVDATA.BinNu); % how many samples spent in a given bin in a given trial
for i = 1:1:LVDATA.BinNu % I need to do it in a complicated way using temporary arrays because if the last trial did not go to the end it gave error (mismatch in TRNu and bin-start at later bins)
    TempArray1 = zeros(LVDATA.TRNu+1,1);
    TempArray1 = find(SampleInBinIsNaN2(:,i)==1)+1; % first sample spend in a given bin, given trial
    if numel(TempArray1)<LVDATA.TRNu+1
        TempArray1(LVDATA.TRNu+1) = NaN;
    end
    EnterIntoBinSampleInd(:,i) = TempArray1(1:LVDATA.TRNu+1);
    TempArray1 = zeros(LVDATA.TRNu+1,1);
    TempArray1 = find(SampleInBinIsNaN3(:,i)==1); % last sample spend in a given bin, given trial
    if numel(TempArray1)<LVDATA.TRNu+1
        TempArray1(LVDATA.TRNu+1) = NaN;
    end
    LeaveBinSampleInd(:,i) = TempArray1(1:LVDATA.TRNu+1); % last sample spend in a given bin, given trial
end
SampleSpentInBin = LeaveBinSampleInd - EnterIntoBinSampleInd + 1; % spent in bin = last sample - first sample +1
LVDATA.EnterIntoBinSampleInd = EnterIntoBinSampleInd;
LVDATA.LeaveBinSampleInd = LeaveBinSampleInd;
LVDATA.SampleSpentInBin = SampleSpentInBin;
% generate longer trial to contain the reward position (which is now in next lap)
HeightPlus10Bin = size(LVDATA.EnterIntoBinSampleInd,1)-2;
LVDATA.EnterIntoBinSampleIndPlusBins = NaN(HeightPlus10Bin,LVDATA.BinNu+10);
LVDATA.EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,1:LVDATA.BinNu) = LVDATA.EnterIntoBinSampleInd(1:HeightPlus10Bin,:);
LVDATA.EnterIntoBinSampleIndPlusBins(1:HeightPlus10Bin,LVDATA.BinNu+1:LVDATA.BinNu+10) = LVDATA.EnterIntoBinSampleInd(2:HeightPlus10Bin+1,1:10);
LVDATA.SampleSpentInBinPlusBins = NaN(HeightPlus10Bin,LVDATA.BinNu+10);
LVDATA.SampleSpentInBinPlusBins(1:HeightPlus10Bin,1:LVDATA.BinNu) = LVDATA.SampleSpentInBin(1:HeightPlus10Bin,:);
LVDATA.SampleSpentInBinPlusBins(1:HeightPlus10Bin,LVDATA.BinNu+1:LVDATA.BinNu+10) = LVDATA.SampleSpentInBin(2:HeightPlus10Bin+1,1:10);

LVDATA.TRNu = sum(LVDATA.EnterIntoBinSampleInd(:,1)>0);

% save behavior plots into BEHAVIOR subfolder
mkdir(sData.sessionInfo.savePath,'Behavior');
savePath = strcat(sData.sessionInfo.savePath,'\Behavior');

% plot HeatBinVelo plot
[LVDATA.VeloBinMatrix, LVDATA.MeanVeloBin]= plotHeatBinVelo(LVDATA,round(LVDATA.stats.VelMaxCmS,-1)); % input= LVDATA.VelDS, Ymax for velo
FileName = strcat('VeloHeatBin-',LVDATA.FileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot HeatBinLicks plot
[LVDATA.LickBinMatrix, LVDATA.MeanLickBin,LVDATA.LickCmMatrix,LVDATA.MeanLickCm] = plotHeatBinLicksPlusBin(LVDATA,LVDATA.BinNu+10); 
FileName = strcat('LickHeatBin-',LVDATA.FileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% plot WaterGiven plot / hit-miss visible
SizeOfRewardZoneCm = 6;
[LVDATA.WaterGivenMatrix,LVDATA.Opto.HitTrialsArray,LVDATA.Opto.HitRate] = plotHeatBinWaterGivenPlusBin(LVDATA,LVDATA.BinNu+10,SizeOfRewardZoneCm);
FileName = strcat('WaterGivenHeatBin-',LVDATA.FileID);
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));



% TRIAL SORTING 
%LVDATA = OptoTrialSorting(LVDATA,TDMSDATA);

%%% Plot stttings
RewardZoneLength = 6 ; % cm
PlusBins = 10;
Xstep = LVDATA.BinSize;
XaxisEnd = (BinNu+PlusBins) * LVDATA.BinSize;
TRNu = LVDATA.TRNu;
LickMax = ceil(max(LVDATA.MeanLickCm));


if IsOpto == 0
    %%% PLOT average lick / trials
    axis = 1:Xstep:XaxisEnd;
    figure('Color','white'); 
    plot(axis,LVDATA.MeanLickCm); hold on
    line([LapLengthSet LapLengthSet],[0 LickMax],'Color','black','LineStyle','--'); hold on   %LickMax
    %line([157+RewardZoneLength 157+RewardZoneLength],[0 LickMax],'Color','black','LineStyle','--');
    xlabel('Position on wheel (cm)');
    ylabel('Licks/cm');
    FileName = strcat('LickSum-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));  
    
    %%% PLOT average of trials speed
    figure('Color','white'); 
    axis = 1:Xstep:BinNu * LVDATA.BinSize;
    plot(axis,LVDATA.MeanVeloBin); 
    xlabel('Position on wheel (cm)');
    ylabel('Velocity (cm/s)');
    FileName = strcat('VeloSum-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
   
elseif IsOpto == 1 

    % TRIAL SORTING 
    LVDATA = OptoTrialSorting2(LVDATA,sData);

    %%% PLOT average lick trial-types
    figure('Color','white'); 
    axis = 1:Xstep:XaxisEnd;
    plot(axis,LVDATA.Opto.MeanLick_LightOff); hold on;
    plot(axis,LVDATA.Opto.MeanLick_LightOn); hold on;
    plot(axis,LVDATA.Opto.MeanLick_LightAfter); hold on;
    line([157 157],[0 LickMax],'Color','black','LineStyle','--');
    hold on
    line([157+RewardZoneLength 157+RewardZoneLength],[0 1.5],'Color','black','LineStyle','--');
    xlabel('Position on wheel (cm)');
    ylabel('Licks/cm');
    legend('Light-Off','Light-On','After-Light','Location','north');
    FileName = strcat('LickSum-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    %%% PLOT Lick Heatmaps
    % LIGHT-OFF
    plotdata = LVDATA.Opto.LickMatrix_LightOff;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    line([157 157],[0 TRNuPlot],'Color','white','LineStyle','--');
    hold on
    line([157+RewardZoneLength 157+RewardZoneLength],[0 TRNu],'Color','white','LineStyle','--');
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-Light-Off'));
    FileName = strcat('LickHeat-LighOff-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % LIGHT-ON
    plotdata = LVDATA.Opto.LickMatrix_LightOn;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    line([157 157],[0 TRNuPlot],'Color','white','LineStyle','--');
    hold on
    line([157+RewardZoneLength 157+RewardZoneLength],[0 TRNu],'Color','white','LineStyle','--');
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-Light-On'));
    FileName = strcat('LickHeat-LighOn-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % AFTER-LIGHT
    plotdata = LVDATA.Opto.LickMatrix_AfterLight;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:XaxisEnd,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 5]);
    j.Label.String = 'Licking (licks/cm)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    line([157 157],[0 TRNuPlot],'Color','white','LineStyle','--');
    hold on
    line([157+RewardZoneLength 157+RewardZoneLength],[0 TRNu],'Color','white','LineStyle','--');
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-AfterLight'));
    FileName = strcat('LickHeat-AfterLigh-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


    %%%% SPEED

    %%% PLOT average of trials speed
    figure('Color','white'); 
    axis = 1:Xstep:BinNu * LVDATA.BinSize;
    plot(axis,LVDATA.Opto.MeanVelo_LightOff); hold on;
    plot(axis,LVDATA.Opto.MeanVelo_LightOn); hold on;
    plot(axis,LVDATA.Opto.MeanVelo_LightAfter); hold on;
    xlabel('Position on wheel (cm)');
    ylabel('Velocity (cm/s)');
    legend('Light-Off','Light-On','After Light','Location','South');
    FileName = strcat('VeloSum-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    %%% PLOT Velo Heatmaps
    % LIGHT-OFF
    plotdata = LVDATA.Opto.VeloMatrix_LightOff;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:BinNu * LVDATA.BinSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-Light-Off'));
    FileName = strcat('VeloHeat-LighOff-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % LIGHT-ON
    plotdata = LVDATA.Opto.VeloMatrix_LightOn;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:BinNu * LVDATA.BinSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-Light-On'));
    FileName = strcat('VeloHeat-LighOn-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

    % AFTER-LIGHT
    plotdata = LVDATA.Opto.VeloMatrix_AfterLight;
    TRNuPlot = size(plotdata,1);
    figure('Color','white'); 
    imagesc(1:Xstep:BinNu * LVDATA.BinSize,1:TRNuPlot,(plotdata)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    j = colorbar; colormap(jet); caxis([0 40]);
    j.Label.String = 'Velocity (cm/s)'; 
    j.Label.FontSize = 11; j.TickDirection = 'out'; 
    xlabel('Position on wheel (cm)'); ylabel('Trials');
    title(strcat(LVDATA.FileID,'-AfterLight'));
    FileName = strcat('VeloHeat-AfterLigh-',LVDATA.FileID);
    savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

end
    
%%%% write data into sData
sData.behavior.wheelPos = AbsDist;                   % double    Absolute position on the wheel
sData.behavior.wheelPosDs = LVDATA.AbsDistDS;        % double    Absolute position on the wheel (downsampled to the imaging FPS)
sData.behavior.wheelPosDsMonIncr = LVDATA.AbsDistDSMonIncr;
%sData.behavior.wheelPosDsBinned                     % double    Absolute position on the wheel (downsampled to the imaging FPS, and binned to some bin size stored in behavior.meta.binSize)
sData.behavior.runSpeed = VelocityInst;              % double    Running speed of the animal
sData.behavior.runSpeedDs = LVDATA.VelDS;            % double    Running speed of the animal (downsampled to the imaging FPS)
sData.behavior.lickDs = LVDATA.LickDS;
sData.behavior.waterRewardDs = LVDATA.WaterRewardDS;
sData.behavior.timeInSec = LVDATA.TimeSec;
sData.behavior.wheelLap = LVDATA.stats.TRNuFull;     % double    Absolute lap number from the start of the recording.
sData.behavior.wheelLapImaging = LVDATA.TRNu;
% Meta
sData.behavior.meta.binSize = LVDATA.BinSize;        % double    Bin size when binning data. For instance, this can be the spatial bin size used for binning linear track data.
sData.behavior.meta.nBins = LVDATA.BinNu;     % double    Number of bins a running wheel has been divided into (useful for place field recordings)
if exist('LVDATA.FakeFrameSignal')
    sData.behavior.meta.fakeframesignal = LVDATA.FakeFrameSignal;
end
sData.behavior.meta.daqSamplingRate = LVDATA.DaqSamplingRate;
sData.behavior.meta.imagingSamplingRate = LVDATA.FrameSignalRate;
% Stats
sData.behavior.stats = LVDATA.stats;
% Details
sData.behavior.details.nFrames = LVDATA.FrameNu;
sData.behavior.details.rewardStartIndices = LVDATA.RewardStartIndex;
sData.behavior.details.frameStartIndices = LVDATA.FrameStartIndex;
sData.behavior.details.hitTrials = LVDATA.Opto.HitTrialsArray;
sData.behavior.details.hitRate = LVDATA.Opto.HitRate;
% Binning
sData.behavior.binning.samplesInBinIndex = LVDATA.SampleInBinMatrix;
sData.behavior.binning.enterIntoBinIndex = LVDATA.EnterIntoBinSampleInd;
sData.behavior.binning.leaveBinIndex = LVDATA.LeaveBinSampleInd;
sData.behavior.binning.samplesSpentInBin = LVDATA.SampleSpentInBin;
sData.behavior.binning.enterIntoBinIndexExtended = LVDATA.EnterIntoBinSampleIndPlusBins;
sData.behavior.binning.samplesSpentInBinExtended = LVDATA.SampleSpentInBinPlusBins;
sData.behavior.binning.veloBinned = LVDATA.VeloBinMatrix;
sData.behavior.binning.lickBinned = LVDATA.LickBinMatrix;
sData.behavior.binning.lickPerCmBinned = LVDATA.LickCmMatrix;
sData.behavior.binning.waterGivenBinned = LVDATA.WaterGivenMatrix;
sData.behavior.binning.meanVeloBinned = LVDATA.MeanVeloBin;
sData.behavior.binning.meanLickBinned = LVDATA.MeanLickBin;
sData.behavior.binning.meanLickPerCmBinned = LVDATA.MeanLickCm;
%Optical stimulation
sData.behavior.opto = LVDATA.Opto;
sData.behavior.opto.IsOptoSession = IsOpto;

% Save file to same path where LV files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end