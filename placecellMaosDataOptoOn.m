function sData = placecellMaosDataOptoOn(sData,type,InOutRatio,ReliabiliyIndex)
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s (set in calcCADATAVelMinNori function input). 
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)
% My data is already filtered for low activity in CaPreProc function(the 99% percentile of ROI data should have been > 3*SD of noise of ROI data)

%%%PARAMATERS TO SET:
%type = 0; %0 use dff data, 1 use deconvolved data, 2 use firing rate data
%IsOpto = 0; % 0 if it is not and optical session, 1 if it is 

MinPlaceFieldSize = 4; % cm minimum size of place field Mao:15
MaxPlaceFieldSize = 120; % cm maximum size of place field Mao:120
%CRITERIA 1 the maximum activity in the potential place field must have been at least 1+activityTreshold -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve
if type == 0 % for dff data
    activityTreshold = 0.6;  
else % deconvolved signal
    activityTreshold = 0.3; % Mao = 0.3
    GaussianFilter = round(sData.behavior.meta.binSize*5); % I set 5 bins to smooth (10 cm) . 2 cm Binsize usually
end
%CRITERIA 2 the mean in-field activity must be at least three times larger than the mean out-of-field activity;
%InOutRatio = 3; % Mao = 3
%CRITERIA 3 more than one third of the trials must have peak position-mapped activity fall within the potential place field
%ReliabiliyIndex = 0.34; % Mao = 0.34

%%% I set place cell tuning around peak to be 60% as larger as minimum, and place field act should be 2x larger than out of field
%%% INITIALIZE PARAMETERS 
% Load parameters:
MaoOpto = struct; % collect place cell data
BinNu = sData.behavior.meta.nBins;
TRNu = sData.behavior.wheelLapImaging;
BinSize = sData.behavior.meta.binSize;
ROINu = sData.imdata.nROIs;


% Create arrays for processed data:
dataset = struct;
if type == 0
    dataset = sData.imdata.binned.RoidFF;
    MaoOpto.datatype = 'dFF';
elseif type == 1
    dataset = sData.imdata.binned.RoiDeconvolved;
    MaoOpto.datatype = 'deconvolved';
elseif type == 2
    dataset = sData.imdata.binned.RoiSpikeRate;
    MaoOpto.datatype = 'spikerate';
end

%{
%%% USE 3-BIN GAUSSIAN SMOOOTHING for deconvolved data. (Mao used 4.5 cm = 3 bin, I will use 3 bin = 6 cm if BinSize is 2 cm)
% Finally not used, since during SpikeRate estimation there is a Gaussian smoothing with FrameRate/2 (15sec), and this new filtering does not cause any effect
if type == 1
    for i = 1:1:ROINu
        DeconvGauss{i} = NaN(TRNu,BinNu); 
        for j = 1:1:TRNu
            SignalTemp1 = dataset{i}(j,:);
            GaussianFiltered = smoothdata(SignalTemp1,'gaussian',5);
            DeconvGauss{i}(j,:) = GaussianFiltered;
        end
    end
    dataset = DeconvGauss;
end
%}

%{
% compare ROI data with and without this Gaussian smoothing:
roi = 29; % set actual roi 
figure();
plot(Mao.PosTuning(roi,:)); hold on;
figure()
plot(GaussianSmoothed(roi,:)); hold on;
caxis([0 inf]);
%}

%%% GENERATE POSITION TUNING CURVES FOR EACH ROIs
% averaging the position-mapped activity across all trials.
LightOnTrials = find(sData.behavior.opto.LightOnTrials==1);
LightOffTrials = find(sData.behavior.opto.LightOffTrials==1);
LightAfterTrials = find(sData.behavior.opto.AfterLightTrials==1);
MaoOpto.LightOff.PosTuning = NaN(ROINu,BinNu);
MaoOpto.LightOn.PosTuning = NaN(ROINu,BinNu);
MaoOpto.LightAfter.PosTuning = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    for j = 1:1:BinNu
        MaoOpto.LightOff.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(LightOffTrials,j));
        MaoOpto.LightOn.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(LightOnTrials,j));  
        MaoOpto.LightAfter.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(LightAfterTrials,j));  
    end
end

%%% USE 3-BIN GAUSSIAN SMOOOTHING for deconvolved data.
if type == 1
    GaussianSmoothedLightOn = NaN(TRNu,BinNu);
    GaussianSmoothedLightOff = NaN(TRNu,BinNu);
    GaussianSmoothedLightAfter = NaN(TRNu,BinNu);
    for i = 1:1:ROINu
        GaussianSmoothedLightOn(i,:) = smoothdata(MaoOpto.LightOn.PosTuningOrig(i,:),'gaussian',GaussianFilter);
        GaussianSmoothedLightOff(i,:) = smoothdata(MaoOpto.LightOff.PosTuningOrig(i,:),'gaussian',GaussianFilter);
        GaussianSmoothedLightAfter(i,:) = smoothdata(MaoOpto.LightAfter.PosTuningOrig(i,:),'gaussian',GaussianFilter);
    end
    MaoOpto.LightOn.PosTuning = GaussianSmoothedLightOn;
    MaoOpto.LightOff.PosTuning = GaussianSmoothedLightOff;
    MaoOpto.LightAfter.PosTuning = GaussianSmoothedLightAfter;
else
    MaoOpto.LightOn.PosTuning = MaoOpto.LightOn.PosTuningOrig;
    MaoOpto.LightOff.PosTuning = MaoOpto.LightOff.PosTuningOrig;
    MaoOpto.LightAfter.PosTuning = MaoOpto.LightAfter.PosTuningOrig;
end


% Normalize mean activity for visualization
MaoOpto.LightOff.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIOff = max(MaoOpto.LightOff.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.LightOff.PosTuningNorm = MaoOpto.LightOff.PosTuning ./ MaxInROIOff; % normalize
MaoOpto.LightOff.PosTuningNorm(isnan(MaoOpto.LightOff.PosTuningNorm)) = 0;

MaoOpto.LightOn.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIOn = max(MaoOpto.LightOn.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.LightOn.PosTuningNorm = MaoOpto.LightOn.PosTuning ./ MaxInROIOn; % normalize
MaoOpto.LightOn.PosTuningNorm(isnan(MaoOpto.LightOn.PosTuningNorm)) = 0;

MaoOpto.LightAfter.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIAfter = max(MaoOpto.LightAfter.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.LightAfter.PosTuningNorm = MaoOpto.LightAfter.PosTuning ./ MaxInROIAfter; % normalize
MaoOpto.LightAfter.PosTuningNorm(isnan(MaoOpto.LightAfter.PosTuningNorm)) = 0;

%{
% PLOT all ROI averaged (all trial activity)
figure('Color','white'); 
imagesc(1:160,1:ROINu,PC_Mao1.SpikeRate_AllTrialsNorm) % For not normalized data plot: ,PC_Mao1.SpikeRate_AllTrials 
c = colorbar;
colormap(jet);
c.Label.String = 'Normalized Spike Rate'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
caxis([0 inf]); %set limits for color plot, below 1st black, above 2nd white
ax = gca; ax.TickDir = 'out'; xticklabels = 0:10:160; xticks = linspace(1, 160, numel(xticklabels));
xlabel('Position on Wheel (cm)');
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)
ylabel('ROIs'); yticklabels = 0:5:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title('Mean activity of all ROIs (all trials, normlized to max)')
fname = 'MeanROIactNorm.fig';
savefig(fullfile(SavePath,fname));
%}

%%% CRITERIA #1: place fields must be a continuous region with MinPlaceFieldSize and MinPlaceFieldSize (cm) in width, within which the activity magnitude 
%  must be above 30% of the difference between the maximum and minimum activity in the position tuning curve;
for k = 1:1:3
    if k == 1
        PosTuning = MaoOpto.LightOff.PosTuning;
    elseif k == 2
        PosTuning = MaoOpto.LightOn.PosTuning;
    elseif k == 3
        PosTuning = MaoOpto.LightAfter.PosTuning;    
    end    
    MinInROI = min(PosTuning,[],2);
    MaxInROI = max(PosTuning,[],2);
    MinInROI(MinInROI<0) = 0;
    Treshold = MinInROI + (MaxInROI - MinInROI)*activityTreshold;   %%% CHANGE 0.3 to 0.6  % activity must be larger in place field than this value (130% of minimum)
    LargerThanTresholdBins = zeros(ROINu,BinNu);
    for i = 1:1:ROINu
        SignalTemp3 = PosTuning(i,:);
        SignalTemp3(SignalTemp3<Treshold(i))=NaN; % put NaN to bins where activity is smaller than treshold
        LargerThanTresholdBins(i,:) = SignalTemp3; % matrix tresholded activity        
    end
    LargerThanTresholdBins2 = ~isnan(LargerThanTresholdBins); % ones when the activity is above treshold. Count ones next to each other
    ExtraZeros = zeros(ROINu,1);
    LargerThanTresholdBins3 = [ExtraZeros LargerThanTresholdBins2 LargerThanTresholdBins2]; % concatenate three matrices to lengten trials 
    LargerThanTresholdBins4 = NaN(ROINu,2*BinNu+1); % calculate potential place field size. Min: 4 cm, max: 100 cm.
    kmax = MaxPlaceFieldSize/BinSize; % max place field is 100 cm, do not check biger size
    for i = 1:1:ROINu
        for j = 1:1:2*BinNu % + kmax is because next trial check
          if LargerThanTresholdBins3(i,j+1) == 1 % find the first / next one, in this matrix data are shifted with one column ->(j+1)
              for m = 1:1:BinNu+kmax-j % count how many ones are after it
                  if LargerThanTresholdBins3(i,j+1+m) == 0 % if there is a zero (low activity), write the length of place field and jump to the next ones           
                      LargerThanTresholdBins4(i,j+1) = m; % how many bin long place field is it , I put an extra column in the beginning
                      break
                  end
              end
          end
       end
    end
    PlaceFieldSize = NaN(ROINu,BinNu);
    for i = 1:1:ROINu
        SignalTemp5 = LargerThanTresholdBins4(i,1:BinNu+1);
        for j = 1:1:BinNu
            if isnan(SignalTemp5(j)) && ~isnan(SignalTemp5(j+1)) 
                PlaceFieldSize(i,j) = SignalTemp5(j+1); % j-1, because it ws sifted in LargerThanTresholdBins matrix to the right
            end
        end
    end
    % Collect ROIs having at least one potential place field
    Criteria1Passed = NaN(ROINu,1);
    FieldMinBin = round(MinPlaceFieldSize/BinSize);
    FieldMaxBin = round(MaxPlaceFieldSize/BinSize);
    for i = 1:1:ROINu
        SignalTemp6 = PlaceFieldSize(i,:);
        if any(SignalTemp6<FieldMaxBin) && any(SignalTemp6>FieldMinBin) % at binsize=2 
            Criteria1Passed(i)=i;
        end
    end

    % How many place fields do they have and in which bin?
    PotPlaceFieldNu = NaN(ROINu,1);
    PotPlaceFieldPos = NaN(ROINu,10);
    PotPlaceFieldLength = NaN(ROINu,10);
    for i = 1:1:ROINu
        counter = 0;
        for j = 1:1:BinNu
            if PlaceFieldSize(i,j) < FieldMaxBin && PlaceFieldSize(i,j) > FieldMinBin % at binsize=2 cm it is 8 and 60 bin
                counter = counter + 1;
                PotPlaceFieldPos(i,counter) = j;
                PotPlaceFieldLength(i,counter) = PlaceFieldSize(i,j);
            end
        end
        PotPlaceFieldNu(i) = counter; 
    end
    if k == 1
        MaoOpto.LightOff.PotPlaceFieldPos = PotPlaceFieldPos;
        MaoOpto.LightOff.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.LightOff.Criteria1Passed = Criteria1Passed;
        MaoOpto.LightOff.Treshold = Treshold;
    elseif k == 2
        MaoOpto.LightOn.PotPlaceFieldPos = PotPlaceFieldPos;
        MaoOpto.LightOn.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.LightOn.Criteria1Passed = Criteria1Passed;
        MaoOpto.LightOn.Treshold = Treshold;
    elseif k == 3
        MaoOpto.LightAfter.PotPlaceFieldPos = PotPlaceFieldPos;
        MaoOpto.LightAfter.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.LightAfter.Criteria1Passed = Criteria1Passed;
        MaoOpto.LightAfter.Treshold = Treshold;
    end
end

%%% CRITERIA 2: the mean in-field activity must be at least three times larger than the mean out-of-field activity;
for k = 1:1:3
    if k == 1  % OFF trials
        PosTuning = MaoOpto.LightOff.PosTuning;
        PotPlaceFieldPos = MaoOpto.LightOff.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightOff.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.LightOff.Criteria1Passed;
    elseif k == 2  % ON trials 
        PosTuning = MaoOpto.LightOn.PosTuning;
        PotPlaceFieldPos = MaoOpto.LightOn.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightOn.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.LightOn.Criteria1Passed;
   elseif k == 3  % ON trials 
        PosTuning = MaoOpto.LightAfter.PosTuning;
        PotPlaceFieldPos = MaoOpto.LightAfter.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightAfter.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.LightAfter.Criteria1Passed;
    end 
    ActAllTrials = [PosTuning PosTuning]; % concatenate two matrices to have to trials consecutive
    ActInPlaceField = NaN(ROINu,10); % max 10 place field might exist
    ActOutPlaceField = NaN(ROINu,10);
    for i = 1:1:ROINu 
        for j = 1:1:10 % if ROI has more than one potential place fields, check all of them
            if isnan(PotPlaceFieldPos(i,j))
                break
            else
            ActInPlaceField(i,j) = mean(ActAllTrials(i,PotPlaceFieldPos(i,j):(PotPlaceFieldPos(i,j)+PotPlaceFieldLength(i,j))));  
            ActOutPlaceField(i,j) = mean(ActAllTrials(i,PotPlaceFieldPos(i,j)+PotPlaceFieldLength(i,j)+1:(PotPlaceFieldPos(i,j)+BinNu)));
            end
        end
    end
    ActRatioInOutPlaceField = ActInPlaceField ./ ActOutPlaceField;
    Criteria12Passed = Criteria1Passed;
    for i = 1:1:ROINu
        if ~any(ActRatioInOutPlaceField(i,:) > InOutRatio) 
           Criteria12Passed(i) = NaN;
        end
    end
    % update potential place field matrix
    if k == 1
        MaoOpto.LightOff.PotPlaceFieldPos2 = MaoOpto.LightOff.PotPlaceFieldPos;
        MaoOpto.LightOff.PotPlaceFieldPos2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightOff.PotPlaceFieldLength2 = MaoOpto.LightOff.PotPlaceFieldLength;
        MaoOpto.LightOff.PotPlaceFieldLength2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightOff.Criteria12Passed = Criteria12Passed;
    elseif k == 2
        MaoOpto.LightOn.PotPlaceFieldPos2 = MaoOpto.LightOn.PotPlaceFieldPos;
        MaoOpto.LightOn.PotPlaceFieldPos2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightOn.PotPlaceFieldLength2 = MaoOpto.LightOn.PotPlaceFieldLength;
        MaoOpto.LightOn.PotPlaceFieldLength2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightOn.Criteria12Passed = Criteria12Passed;
    elseif k == 3
        MaoOpto.LightAfter.PotPlaceFieldPos2 = MaoOpto.LightAfter.PotPlaceFieldPos;
        MaoOpto.LightAfter.PotPlaceFieldPos2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightAfter.PotPlaceFieldLength2 = MaoOpto.LightAfter.PotPlaceFieldLength;
        MaoOpto.LightAfter.PotPlaceFieldLength2(ActRatioInOutPlaceField<3) = NaN;
        MaoOpto.LightAfter.Criteria12Passed = Criteria12Passed;
    end
end

%%% CRITERIA 3: more than one third of the trials must have peak position-mapped activity fall within the potential place field
for k = 1:1:3  
    Reliability = NaN(ROINu,1);
    if k == 1  % OFF trials
        PotPlaceFieldPos2 = MaoOpto.LightOff.PotPlaceFieldPos2;
        PotPlaceFieldLength2 = MaoOpto.LightOff.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.LightOff.Criteria12Passed;
        TrialNu = length(LightOffTrials);
        Trials = LightOffTrials;
    elseif k == 2
        PotPlaceFieldPos2 = MaoOpto.LightOn.PotPlaceFieldPos2;
        PotPlaceFieldLength2 = MaoOpto.LightOn.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.LightOn.Criteria12Passed;
        TrialNu = length(LightOnTrials);
        Trials = LightOnTrials;
    elseif k == 3
        PotPlaceFieldPos2 = MaoOpto.LightAfter.PotPlaceFieldPos2;
        PotPlaceFieldLength2 = MaoOpto.LightAfter.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.LightAfter.Criteria12Passed;
        TrialNu = length(LightAfterTrials);
        Trials = LightAfterTrials;
    end
    for i = 1:1:ROINu
        if any(~isnan(PotPlaceFieldPos2(i,:))) % if no potential place field, jump to next roi
            counter = 0; % counter for reliable trials
            PFStart = PotPlaceFieldPos2(i); % place field start (bin)
            PFEnd = PotPlaceFieldPos2(i) + PotPlaceFieldLength2(i); % place field end
            for m = 1:1:TrialNu
                TrialData = dataset{1,i}(Trials(m),:); % put each trial data into a container
                if any(isnan(TrialData)) % after the end  of recording there are NaNs. If this is the end, put reliability data into array and Go to next ROI
                   break 
                end
                if mean(TrialData) == 0 % sometimes there are only zeros in a row. it makes matlab crazy in the next trial, jump ot next trial
                    continue
                end
                MaxBin = find(TrialData == max(TrialData(:))); % find maximum (which bin)
                MaxBinFirst = MaxBin(1,1); % sometimes the same max data is in two neighboring position, take the first
                if PFEnd <= BinNu
                    if MaxBinFirst >= PFStart && MaxBinFirst <= PFEnd % peak should be between place field start and end. Place field might contimu in to next bin
                        counter = counter + 1;
                    end
                elseif PFEnd > BinNu % if place field goes  through trial start
                    if MaxBinFirst >= PFStart || MaxBinFirst <= rem(PFEnd,BinNu)
                        counter = counter + 1;
                    end
                end
            end
            Reliability(i,1) = counter / TrialNu;  % reliability data, 1 is the max. Needed >0.33 to be a place cell
        end
    end
    Criteria123Passed = Criteria12Passed;
    for i = 1:1:ROINu
        if Reliability(i,1) < ReliabiliyIndex
            Criteria123Passed(i) = NaN;
        end
    end
    PlaceCells =  Criteria123Passed(~isnan(Criteria123Passed));
    % update potential place field matrix 
    PotPlaceFieldPos3 = PotPlaceFieldPos2;
    PotPlaceFieldPos3(Reliability < ReliabiliyIndex) = NaN;
    PotPlaceFieldLength3 = PotPlaceFieldLength2;
    PotPlaceFieldLength3(Reliability < ReliabiliyIndex) = NaN;
    if k == 1
        MaoOpto.LightOff.PlaceFieldStartBin = PotPlaceFieldPos3;
        MaoOpto.LightOff.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.LightOff.Criteria123Passed = Criteria123Passed;
        MaoOpto.LightOff.Reliability = Reliability;
        MaoOpto.LightOff.PlaceCells = PlaceCells;
    elseif k == 2
        MaoOpto.LightOn.PlaceFieldStartBin = PotPlaceFieldPos3;
        MaoOpto.LightOn.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.LightOn.Criteria123Passed = Criteria123Passed;
        MaoOpto.LightOn.Reliability = Reliability;
        MaoOpto.LightOn.PlaceCells = PlaceCells;
    elseif k == 3
        MaoOpto.LightAfter.PlaceFieldStartBin = PotPlaceFieldPos3;
        MaoOpto.LightAfter.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.LightAfter.Criteria123Passed = Criteria123Passed;
        MaoOpto.LightAfter.Reliability = Reliability;
        MaoOpto.LightAfter.PlaceCells = PlaceCells;
    end
end


% PLOT place cell ROI averaged (all trial activity)
MaoOpto.LightOff.placeROINu = length(MaoOpto.LightOff.PlaceCells);  
MaoOpto.LightOff.placeROINormActBin = MaoOpto.LightOff.PosTuningNorm(MaoOpto.LightOff.PlaceCells,:); % collect only place cell data
MaoOpto.LightOn.placeROINu = length(MaoOpto.LightOn.PlaceCells);  
MaoOpto.LightOn.placeROINormActBin = MaoOpto.LightOn.PosTuningNorm(MaoOpto.LightOn.PlaceCells,:); % collect only place cell data
MaoOpto.LightAfter.placeROINu = length(MaoOpto.LightAfter.PlaceCells);  
MaoOpto.LightAfter.placeROINormActBin = MaoOpto.LightAfter.PosTuningNorm(MaoOpto.LightAfter.PlaceCells,:); % collect only place cell data

if type == 0
    sData.imdata.MaoPC_dff = struct;
    sData.imdata.MaoPC_dff = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell-dff');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell-dff');
elseif type == 1
    sData.imdata.MaoPC_deconv = struct;
    sData.imdata.MaoPC_deconv = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell-deconv');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell-deconv');
elseif type == 2
    sData.imdata.MaoPC_deconv = struct;
    sData.imdata.MaoPC_deconv = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell-spikerate');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell-spikerate');
end
% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% plot position tuning curves with potential place fields
% plot only place cells:
for k = 1:1:3
    if k == 1 
        PlaceCells = MaoOpto.LightOff.PlaceCells;
        PosTuning = MaoOpto.LightOff.PosTuning;
        Treshold = MaoOpto.LightOff.Treshold;
        PlaceFieldStartBin = MaoOpto.LightOff.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.LightOff.PlaceFieldBinLength;
        PotPlaceFieldPos = MaoOpto.LightOff.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightOff.PotPlaceFieldLength;
        text = 'Laser-off';
        file = '-postun-LaserOff';
        mkdir(savePathPre,'\MaoPlaceCellsLaserOff');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsLaserOff');
    elseif k == 2
        PlaceCells = MaoOpto.LightOn.PlaceCells;
        PosTuning = MaoOpto.LightOn.PosTuning;
        Treshold = MaoOpto.LightOn.Treshold;
        PlaceFieldStartBin = MaoOpto.LightOn.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.LightOn.PlaceFieldBinLength;
        PotPlaceFieldPos = MaoOpto.LightOn.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightOn.PotPlaceFieldLength;
        text = 'Laser-on';
        file = '-postun-LaserOn';
        mkdir(savePathPre,'\MaoPlaceCellsLaserOn');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsLaserOn');
    elseif k == 3
        PlaceCells = MaoOpto.LightAfter.PlaceCells;
        PosTuning = MaoOpto.LightAfter.PosTuning;
        Treshold = MaoOpto.LightAfter.Treshold;
        PlaceFieldStartBin = MaoOpto.LightAfter.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.LightAfter.PlaceFieldBinLength;
        PotPlaceFieldPos = MaoOpto.LightAfter.PotPlaceFieldPos;
        PotPlaceFieldLength = MaoOpto.LightAfter.PotPlaceFieldLength;
        text = 'After-laser';
        file = '-postun-deconv-LaserAfter';
        mkdir(savePathPre,'\MaoPlaceCellsAfterLaser');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsAfterLaser');
    end
    for j = 1:1:length(PlaceCells)
        i = PlaceCells(j);
        figure();
        Xaxis = BinSize:BinSize:BinSize*BinNu;
        Ymax = (max(PosTuning(i,:)))+0.00000000001;
        rectangle('Position',[PlaceFieldStartBin(i,1)*BinSize Ymax/500 PlaceFieldBinLength(i,1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
        plot(Xaxis,PosTuning(i,1:BinNu),'LineWidth',2); hold on; % pos tuning curve
        line([0 BinSize*BinNu],[Treshold(i) Treshold(i)],'Color','red','LineWidth',1); hold on; 
        line([PotPlaceFieldPos(i,1)*BinSize PotPlaceFieldPos(i,1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; line([(PotPlaceFieldPos(i,1)+PotPlaceFieldLength(i,1))*BinSize (PotPlaceFieldPos(i,1)+PotPlaceFieldLength(i,1))*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
        line([PotPlaceFieldPos(i,2)*BinSize PotPlaceFieldPos(i,2)*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on; line([(PotPlaceFieldPos(i,2)+PotPlaceFieldLength(i,2))*BinSize (PotPlaceFieldPos(i,2)+PotPlaceFieldLength(i,2))*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on;
        axis([0 160 0 Ymax]); % ceil(Ymax)
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(i),'-',text));
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning of activity');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(i),file);
        savefig(fullfile(SavePath,fname));
        saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
   end
   close all
end


% sort place cell activity based on peak activity
GaussianFilter = 5; 
plotSortROIsMaxActHeatsDataPlaceOpto(sData,GaussianFilter,type);

end
