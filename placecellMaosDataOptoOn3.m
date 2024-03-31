function sData = placecellMaosDataOptoOn3(sData,datatype,FigVisible)
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)

%CRITERIA 1 the maximum activity in the potential place field must have been at least 1+activityTreshold -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve
%CRITERIA 2 the mean in-field activity must be at least three times larger than the mean out-of-field activity;
%CRITERIA 3 more than one third of the trials must have peak position-mapped activity fall within the potential place field


%%%PARAMATERS TO SET:
%type = 0; %0 use dff data, 1 use deconvolved data, 2 use firing rate data
%IsOpto = 0; % 0 if it is not and optical session, 1 if it is 

MinPlaceFieldSize = 2; % cm minimum size of place field Mao:15
MaxPlaceFieldSize = 100; % cm maximum size of place field Mao:120
%CRITERIA 1 the maximum activity in the potential place field must have been at least 1+activityTreshold -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve
if datatype == 0 % for dff data
    activityTreshold = 0.5;  
    InOutRatio = 2;
    ReliabiliyIndex = 0.3;
else % deconvolved signal
    activityTreshold = 0.3; % Mao = 0.3
    InOutRatio = 3; % Mao = 3
    ReliabiliyIndex = 0.2; % Mao = 0.34
    GaussianFilter = round(sData.behavior.meta.binSize*5); % I set 5 bins to smooth (10 cm) . 2 cm Binsize usually
end

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
if datatype == 0
    dataset = sData.imdata.binned.RoidFF;
    MaoOpto.datatype = 'dFF';
elseif datatype == 1
    dataset = sData.imdata.binned.RoiDeconvolved;
    MaoOpto.datatype = 'deconvolved';
elseif datatype == 2
    dataset = sData.imdata.binned.RoiSpikeRate;
    MaoOpto.datatype = 'spikerate';
end

%%% GENERATE POSITION TUNING CURVES FOR EACH ROIs
% averaging the position-mapped activity across all trials.
OptoOnTrials = sData.imdata.binned.OptoOnTrials;
OptoOffTrials = sData.imdata.binned.OptoOffTrials;
OptoAfterTrials = sData.imdata.binned.AfterOptoTrials;
MaoOpto.OptoOff.PosTuning = NaN(ROINu,BinNu);
MaoOpto.OptoOn.PosTuning = NaN(ROINu,BinNu);
MaoOpto.OptoAfter.PosTuning = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    for j = 1:1:BinNu
        MaoOpto.OptoOff.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(OptoOffTrials,j));
        MaoOpto.OptoOn.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(OptoOnTrials,j));  
        MaoOpto.OptoAfter.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(OptoAfterTrials,j));  
    end
end

%%% USE 3-BIN GAUSSIAN SMOOOTHING for deconvolved data.
if datatype == 1
    GaussianSmoothedOptoOn = NaN(TRNu,BinNu);
    GaussianSmoothedOptoOff = NaN(TRNu,BinNu);
    GaussianSmoothedOptoAfter = NaN(TRNu,BinNu);
    for i = 1:1:ROINu
        GaussianSmoothedOptoOn(i,:) = smoothdata(MaoOpto.OptoOn.PosTuningOrig(i,:),'gaussian',GaussianFilter);
        GaussianSmoothedOptoOff(i,:) = smoothdata(MaoOpto.OptoOff.PosTuningOrig(i,:),'gaussian',GaussianFilter);
        GaussianSmoothedOptoAfter(i,:) = smoothdata(MaoOpto.OptoAfter.PosTuningOrig(i,:),'gaussian',GaussianFilter);
    end
    MaoOpto.OptoOn.PosTuning = GaussianSmoothedOptoOn;
    MaoOpto.OptoOff.PosTuning = GaussianSmoothedOptoOff;
    MaoOpto.OptoAfter.PosTuning = GaussianSmoothedOptoAfter;
else
    MaoOpto.OptoOn.PosTuning = MaoOpto.OptoOn.PosTuningOrig;
    MaoOpto.OptoOff.PosTuning = MaoOpto.OptoOff.PosTuningOrig;
    MaoOpto.OptoAfter.PosTuning = MaoOpto.OptoAfter.PosTuningOrig;
end


% Normalize mean activity for visualization
MaoOpto.OptoOff.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIOff = max(MaoOpto.OptoOff.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.OptoOff.PosTuningNorm = MaoOpto.OptoOff.PosTuning ./ MaxInROIOff; % normalize
MaoOpto.OptoOff.PosTuningNorm(isnan(MaoOpto.OptoOff.PosTuningNorm)) = 0;

MaoOpto.OptoOn.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIOn = max(MaoOpto.OptoOn.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.OptoOn.PosTuningNorm = MaoOpto.OptoOn.PosTuning ./ MaxInROIOn; % normalize
MaoOpto.OptoOn.PosTuningNorm(isnan(MaoOpto.OptoOn.PosTuningNorm)) = 0;

MaoOpto.OptoAfter.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIAfter = max(MaoOpto.OptoAfter.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.OptoAfter.PosTuningNorm = MaoOpto.OptoAfter.PosTuning ./ MaxInROIAfter; % normalize
MaoOpto.OptoAfter.PosTuningNorm(isnan(MaoOpto.OptoAfter.PosTuningNorm)) = 0;


%%% CRITERIA #1: place fields must be a continuous region with MinPlaceFieldSize and MinPlaceFieldSize (cm) in width, within which the activity magnitude 
%  must be above X% of the difference between the maximum and minimum activity in the position tuning curve;
for k = 1:1:3
    if k == 1
        PosTuning = MaoOpto.OptoOff.PosTuning;
    elseif k == 2
        PosTuning = MaoOpto.OptoOn.PosTuning;
    elseif k == 3
        PosTuning = MaoOpto.OptoAfter.PosTuning;    
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
        if ~isnan(SignalTemp5(2)) && ~isnan(SignalTemp5(BinNu+1)) % if there is activity at the end and beginning of track as well
            % delete values from the beginning of the track if there is already high activity at the end
            SignalTemp5(2:SignalTemp5(2)+1) = 0;
        end
        for j = 1:1:BinNu
            if isnan(SignalTemp5(j)) && ~isnan(SignalTemp5(j+1)) 
                PlaceFieldSize(i,j) = SignalTemp5(j+1); % j+1, because it was sifted in LargerThanTresholdBins matrix to the right
            end
        end
    end
    % Collect ROIs having at least one potential place field
    Criteria1Passed = NaN(ROINu,1);
    FieldMinBin = round(MinPlaceFieldSize/BinSize);
    FieldMaxBin = round(MaxPlaceFieldSize/BinSize);
    for i = 1:1:ROINu
        SignalTemp6 = PlaceFieldSize(i,:);
        if any(SignalTemp6 < FieldMaxBin) && any(SignalTemp6 > FieldMinBin)  
            Criteria1Passed(i)=i;
        end
    end

    % How many place fields do they have and in which bin?
    PotPlaceFieldNu = NaN(ROINu,1);
    PotPlaceFieldStartBin = NaN(ROINu,20);
    PotPlaceFieldLength = NaN(ROINu,20);
    for i = 1:1:ROINu
        counter = 0;
        for j = 1:1:BinNu
            if PlaceFieldSize(i,j) < FieldMaxBin && PlaceFieldSize(i,j) > FieldMinBin % at binsize=2 cm it is 8 and 60 bin
                counter = counter + 1;
                PotPlaceFieldStartBin(i,counter) = j;
                PotPlaceFieldLength(i,counter) = PlaceFieldSize(i,j);
            end
        end
        PotPlaceFieldNu(i) = counter; 
    end
    if k == 1
        MaoOpto.OptoOff.PotPlaceFieldStartBin = PotPlaceFieldStartBin;
        MaoOpto.OptoOff.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.OptoOff.Criteria1Passed = Criteria1Passed;
        MaoOpto.OptoOff.Treshold = Treshold;
    elseif k == 2
        MaoOpto.OptoOn.PotPlaceFieldStartBin = PotPlaceFieldStartBin;
        MaoOpto.OptoOn.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.OptoOn.Criteria1Passed = Criteria1Passed;
        MaoOpto.OptoOn.Treshold = Treshold;
    elseif k == 3
        MaoOpto.OptoAfter.PotPlaceFieldStartBin = PotPlaceFieldStartBin;
        MaoOpto.OptoAfter.PotPlaceFieldLength = PotPlaceFieldLength;
        MaoOpto.OptoAfter.Criteria1Passed = Criteria1Passed;
        MaoOpto.OptoAfter.Treshold = Treshold;
    end
end

%%% CRITERIA 2: the mean in-field activity must be at least X-times larger than the mean out-of-field activity;
for k = 1:1:3
    if k == 1  % OFF trials
        PosTuning = MaoOpto.OptoOff.PosTuning;
        PotPlaceFieldStartBin = MaoOpto.OptoOff.PotPlaceFieldStartBin;
        PotPlaceFieldLength = MaoOpto.OptoOff.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.OptoOff.Criteria1Passed;
    elseif k == 2  % ON trials 
        PosTuning = MaoOpto.OptoOn.PosTuning;
        PotPlaceFieldStartBin = MaoOpto.OptoOn.PotPlaceFieldStartBin;
        PotPlaceFieldLength = MaoOpto.OptoOn.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.OptoOn.Criteria1Passed;
   elseif k == 3  % ON trials 
        PosTuning = MaoOpto.OptoAfter.PosTuning;
        PotPlaceFieldStartBin = MaoOpto.OptoAfter.PotPlaceFieldStartBin;
        PotPlaceFieldLength = MaoOpto.OptoAfter.PotPlaceFieldLength;
        Criteria1Passed = MaoOpto.OptoAfter.Criteria1Passed;
    end 
    ActAllTrials = [PosTuning PosTuning]; % concatenate two matrices to have to trials consecutive
    ActInPlaceField = NaN(ROINu,10); % max 10 place field might exist
    ActOutPlaceField = NaN(ROINu,10);
    for i = 1:1:ROINu 
        for j = 1:1:20 % if ROI has more than one potential place fields, check all of them
            if isnan(PotPlaceFieldStartBin(i,j))
                break
            else
            ActInPlaceField(i,j) = mean(ActAllTrials(i,PotPlaceFieldStartBin(i,j):(PotPlaceFieldStartBin(i,j)+PotPlaceFieldLength(i,j)-1)));  
            ActOutPlaceField(i,j) = mean(ActAllTrials(i,PotPlaceFieldStartBin(i,j)+PotPlaceFieldLength(i,j):(PotPlaceFieldStartBin(i,j)+BinNu)-1));
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
        MaoOpto.OptoOff.PotPlaceFieldStartBin2 = MaoOpto.OptoOff.PotPlaceFieldStartBin;
        MaoOpto.OptoOff.PotPlaceFieldStartBin2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoOff.PotPlaceFieldLength2 = MaoOpto.OptoOff.PotPlaceFieldLength;
        MaoOpto.OptoOff.PotPlaceFieldLength2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoOff.Criteria12Passed = Criteria12Passed;
    elseif k == 2
        MaoOpto.OptoOn.PotPlaceFieldStartBin2 = MaoOpto.OptoOn.PotPlaceFieldStartBin;
        MaoOpto.OptoOn.PotPlaceFieldStartBin2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoOn.PotPlaceFieldLength2 = MaoOpto.OptoOn.PotPlaceFieldLength;
        MaoOpto.OptoOn.PotPlaceFieldLength2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoOn.Criteria12Passed = Criteria12Passed;
    elseif k == 3
        MaoOpto.OptoAfter.PotPlaceFieldStartBin2 = MaoOpto.OptoAfter.PotPlaceFieldStartBin;
        MaoOpto.OptoAfter.PotPlaceFieldStartBin2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoAfter.PotPlaceFieldLength2 = MaoOpto.OptoAfter.PotPlaceFieldLength;
        MaoOpto.OptoAfter.PotPlaceFieldLength2(ActRatioInOutPlaceField < InOutRatio) = NaN;
        MaoOpto.OptoAfter.Criteria12Passed = Criteria12Passed;
    end
end

%%% CRITERIA 3: more than one third of the trials must have peak position-mapped activity fall within the potential place field
for k = 1:1:3  
    Reliability = NaN(ROINu,1);
    if k == 1  % OFF trials
        PotPlaceFieldStartBin2 = MaoOpto.OptoOff.PotPlaceFieldStartBin2;
        PotPlaceFieldLength2 = MaoOpto.OptoOff.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.OptoOff.Criteria12Passed;
        TrialNu = length(OptoOffTrials);
        Trials = OptoOffTrials;
    elseif k == 2
        PotPlaceFieldStartBin2 = MaoOpto.OptoOn.PotPlaceFieldStartBin2;
        PotPlaceFieldLength2 = MaoOpto.OptoOn.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.OptoOn.Criteria12Passed;
        TrialNu = length(OptoOnTrials);
        Trials = OptoOnTrials;
    elseif k == 3
        PotPlaceFieldStartBin2 = MaoOpto.OptoAfter.PotPlaceFieldStartBin2;
        PotPlaceFieldLength2 = MaoOpto.OptoAfter.PotPlaceFieldLength2;
        Criteria12Passed = MaoOpto.OptoAfter.Criteria12Passed;
        TrialNu = length(OptoAfterTrials);
        Trials = OptoAfterTrials;
    end
    for i = 1:1:ROINu
        if any(~isnan(PotPlaceFieldStartBin2(i,:))) % if no potential place field, jump to next roi
            counter = 0; % counter for reliable trials
            %%% 
            PosTuningPeakSearch = NaN(1,BinNu); % in position tuning curve mask out those parts where no acceptable peak was found  
            PotBins = PotPlaceFieldStartBin2(i,:); % look for potential PF positions after CRITERIA 2 is tested
            PotBins = PotBins(~isnan(PotBins));
            PotLength = PotPlaceFieldLength2(i,:);
            PotLength = PotLength(~isnan(PotLength));
            PotPlaceFieldNu = sum(~isnan(PotPlaceFieldStartBin2(i,:)));
            for j = 1:1:PotPlaceFieldNu
                PosTuningPeakSearch(PotBins(j):(PotBins(j)+PotLength(j)-1)) = ActAllTrials(i,(PotBins(j):PotBins(j)+PotLength(j)-1)); % get data from circularized Position Tuning dataset
            end
            PeakMax = max(PosTuningPeakSearch); % if more peaks pass the criteria, choose the larger peak
            PeakMaxLoc = find(PosTuningPeakSearch==PeakMax,1); % find the position (bin) 
            PotBins(PotBins > PeakMaxLoc) = NaN;
            PotBins = PotBins(~isnan(PotBins));
            PFStart = max(PotBins); % search to which potential place field the max peak belongs to 
            PFLength = PotPlaceFieldLength2(i,PotPlaceFieldStartBin2(i,:)==PFStart);
            PFEnd = PFStart + PFLength-1; % place field end
            %%%
            %PFStart = PotPlaceFieldStartBin2(i); % place field start (bin)
            %PFEnd = PotPlaceFieldStartBin2(i) + PotPlaceFieldLength2(i)-1; % place field end           
            for m = 1:1:TrialNu
                TrialData = dataset{1,i}(Trials(m),:); % put each trial data into a container
                if any(isnan(TrialData)) % after the end  of recording there are NaNs. If this is the end, put reliability data into array and Go to next ROI
                   break 
                end
                if mean(TrialData) == 0 % sometimes there are only zeros in a row. it makes matlab crazy in the next trial, jump ot next trial
                    continue
                end
                MaxBinFirst = find(TrialData == max(TrialData(:)),1); % find maximum (which bin) sometimes the same max data is in two neighboring position, take the first
                MaxSmoothed = max(smoothdata(TrialData(:),'movmean',3));
                if MaxSmoothed < (mean(TrialData) + 2*(std(TrialData))) % the peak has to be larger than to mean+noise to be consideres as a peak, otherwise it is considered as noise
                    continue
                end                
                
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
    PotPlaceFieldStartBin3 = PotPlaceFieldStartBin2;
    PotPlaceFieldStartBin3(Reliability < ReliabiliyIndex) = NaN;
    PotPlaceFieldLength3 = PotPlaceFieldLength2;
    PotPlaceFieldLength3(Reliability < ReliabiliyIndex) = NaN;
    if k == 1
        MaoOpto.OptoOff.PlaceFieldStartBin = PotPlaceFieldStartBin3;
        MaoOpto.OptoOff.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.OptoOff.Criteria123Passed = Criteria123Passed;
        MaoOpto.OptoOff.Reliability = Reliability;
        MaoOpto.OptoOff.PlaceCells = PlaceCells;
    elseif k == 2
        MaoOpto.OptoOn.PlaceFieldStartBin = PotPlaceFieldStartBin3;
        MaoOpto.OptoOn.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.OptoOn.Criteria123Passed = Criteria123Passed;
        MaoOpto.OptoOn.Reliability = Reliability;
        MaoOpto.OptoOn.PlaceCells = PlaceCells;
    elseif k == 3
        MaoOpto.OptoAfter.PlaceFieldStartBin = PotPlaceFieldStartBin3;
        MaoOpto.OptoAfter.PlaceFieldBinLength = PotPlaceFieldLength3;
        MaoOpto.OptoAfter.Criteria123Passed = Criteria123Passed;
        MaoOpto.OptoAfter.Reliability = Reliability;
        MaoOpto.OptoAfter.PlaceCells = PlaceCells;
    end
end


% PLOT place cell ROI averaged (all trial activity)
MaoOpto.OptoOff.placeROINu = length(MaoOpto.OptoOff.PlaceCells);  
MaoOpto.OptoOff.placeROINormActBin = MaoOpto.OptoOff.PosTuningNorm(MaoOpto.OptoOff.PlaceCells,:); % collect only place cell data
MaoOpto.OptoOn.placeROINu = length(MaoOpto.OptoOn.PlaceCells);  
MaoOpto.OptoOn.placeROINormActBin = MaoOpto.OptoOn.PosTuningNorm(MaoOpto.OptoOn.PlaceCells,:); % collect only place cell data
MaoOpto.OptoAfter.placeROINu = length(MaoOpto.OptoAfter.PlaceCells);  
MaoOpto.OptoAfter.placeROINormActBin = MaoOpto.OptoAfter.PosTuningNorm(MaoOpto.OptoAfter.PlaceCells,:); % collect only place cell data

MaoOpto.PlaceCellinAllProt = intersect(intersect(MaoOpto.OptoOff.PlaceCells,MaoOpto.OptoOn.PlaceCells),MaoOpto.OptoAfter.PlaceCells);%,MaoOpto.OptoAfter.PlaceCells);
MaoOpto.PlaceCellinAnyProt = union(union(MaoOpto.OptoOff.PlaceCells,MaoOpto.OptoOn.PlaceCells),MaoOpto.OptoAfter.PlaceCells);%,MaoOpto.OptoAfter.PlaceCells);

if datatype == 0
    sData.imdata.MaoPC_Opto_dff_new = struct;
    sData.imdata.MaoPC_Opto_dff_new = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCellOpto-dff-new');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCellOpto-dff-new');
elseif datatype == 1
    sData.imdata.MaoPC_Opto_deconv_new = struct;
    sData.imdata.MaoPC_Opto_deconv_new = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCellOpto-deconv-new');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCellOpto-deconv-new');
elseif datatype == 2
    sData.imdata.MaoPC_Opto_spikerate_new = struct;
    sData.imdata.MaoPC_Opto_spikerate_new = MaoOpto;
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCellOpto-spikerate-new');
    savePathPre = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCellOpto-spikerate-new');
end

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% plot position tuning curves with potential place fields
% plot only place cells:
for k = 1:1:3
    if k == 1 
        PlaceCells = MaoOpto.OptoOff.PlaceCells;
        PosTuning = MaoOpto.OptoOff.PosTuning;
        Treshold = MaoOpto.OptoOff.Treshold;
        PlaceFieldStartBin = MaoOpto.OptoOff.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.OptoOff.PlaceFieldBinLength;
        text = 'Opto-off';
        file = '-postun-OptoOff';
        mkdir(savePathPre,'\MaoPlaceCellsOptoOff');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsOptoOff');
    elseif k == 2
        PlaceCells = MaoOpto.OptoOn.PlaceCells;
        PosTuning = MaoOpto.OptoOn.PosTuning;
        Treshold = MaoOpto.OptoOn.Treshold;
        PlaceFieldStartBin = MaoOpto.OptoOn.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.OptoOn.PlaceFieldBinLength;
        text = 'Opto-on';
        file = '-postun-OptoOn';
        mkdir(savePathPre,'\MaoPlaceCellsOptoOn');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsOptoOn');
    elseif k == 3
        PlaceCells = MaoOpto.OptoAfter.PlaceCells;
        PosTuning = MaoOpto.OptoAfter.PosTuning;
        Treshold = MaoOpto.OptoAfter.Treshold;
        PlaceFieldStartBin = MaoOpto.OptoAfter.PlaceFieldStartBin;
        PlaceFieldBinLength = MaoOpto.OptoAfter.PlaceFieldBinLength;
        text = 'After-Opto';
        file = '-postun-deconv-OptoAfter';
        mkdir(savePathPre,'\MaoPlaceCellsAfterOpto');
        SavePath = strcat(savePathPre,'\MaoPlaceCellsAfterOpto');
    end
    for j = 1:1:length(PlaceCells)
        i = PlaceCells(j);
        figure('visible',FigVisible);
        Xaxis = BinSize:BinSize:BinSize*BinNu;
        Ymax = (max(PosTuning(i,:)))+0.00000000001;
        for k = 1:1:size(PlaceFieldStartBin,2)
            if ~isnan(PlaceFieldStartBin(i,k))
                if PlaceFieldStartBin(i,k) + PlaceFieldBinLength(i,k)-1 < BinNu +1
                    rectangle('Position',[PlaceFieldStartBin(i,k)*BinSize Ymax/500 (PlaceFieldBinLength(i,k)-1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
                    line([PlaceFieldStartBin(i,k)*BinSize PlaceFieldStartBin(i,k)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; line([(PlaceFieldStartBin(i,k)+PlaceFieldBinLength(i,k)-1)*BinSize (PlaceFieldStartBin(i,k)+PlaceFieldBinLength(i,k)-1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
                else
                    rectangle('Position',[PlaceFieldStartBin(i,k)*BinSize Ymax/500 (BinNu-PlaceFieldStartBin(i,k))*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
                    rectangle('Position',[BinSize Ymax/500 (PlaceFieldStartBin(i,k)+PlaceFieldBinLength(i,k)-BinNu-1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
                    line([PlaceFieldStartBin(i,k)*BinSize PlaceFieldStartBin(i,k)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; 
                    line([(PlaceFieldStartBin(i,k)+PlaceFieldBinLength(i,k)-BinNu-1)*BinSize (PlaceFieldStartBin(i,k)+PlaceFieldBinLength(i,k)-BinNu-1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
                end
            end
        end
        plot(Xaxis,PosTuning(i,1:BinNu),'LineWidth',2); hold on; % pos tuning curve
        line([0 BinSize*BinNu],[Treshold(i) Treshold(i)],'Color','red','LineWidth',1); hold on; 
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
plotSortROIsMaxActHeatsDataPlaceOpto2(sData,GaussianFilter,datatype,savePathPre);

end
