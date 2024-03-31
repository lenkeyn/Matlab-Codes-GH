function sData = placecellMaosDataOptoOnMoreProt(sData,type,InOutRatio,ReliabiliyIndex,Protocol)
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s (set in calcCADATAVelMinNori function input). 
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)
% My data is already filtered for low activity in CaPreProc function(the 99% percentile of ROI data should have been > 3*SD of noise of ROI data)

%%%PARAMATERS TO SET:
%type = 0; %0 use dff data, 1 use deconvolved data, 2 use firing rate data
%IsOpto = 0; % 0 if it is not and optical session, 1 if it is 
%GaussianFilter = 5; 

MinPlaceFieldSize = 2; % cm minimum size of place field Mao:15, I used 4 cm for analysis until 2021.09., then 2 cm
MaxPlaceFieldSize = 100; % cm maximum size of place field Mao:120, I used 80 cm for PCs for analysis until 2021.09., then 100
%CRITERIA 1 the maximum activity in the potential place field must have been at least 1+activityTreshold -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve
if type == 0 % for dff data
    activityTreshold = 0.5;   % I used 0.6 for analysis until 2021.09., then 0.5
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

%%% GENERATE POSITION TUNING CURVES FOR EACH ROIs
% averaging the position-mapped activity across all trials.
if exist  sData.behavior.optoMoreProts.OptoStimProtTrials
Trials = find(sData.behavior.optoMoreProts.OptoStimProtTrials==Protocol);
if exist  sData.behavior.optoMoreProts.OptoStimProtTrialsReal
MaoOpto.PosTuning = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    for j = 1:1:BinNu
        MaoOpto.PosTuningOrig(i,j) =  nanmean(dataset{1,i}(Trials,j));  
    end
end

%%% USE 3-BIN GAUSSIAN SMOOOTHING for deconvolved data.
if type == 1
    GaussianSmoothed = NaN(TRNu,BinNu);
    for i = 1:1:ROINu
        GaussianSmoothed(i,:) = smoothdata(MaoOpto.PosTuningOrig(i,:),'gaussian',GaussianFilter);
    end
    MaoOpto.PosTuning = GaussianSmoothed;
else
    MaoOpto.PosTuning = MaoOpto.PosTuningOrig;
end


% Normalize mean activity for visualization
MaoOpto.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROIOn = max(MaoOpto.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
MaoOpto.PosTuningNorm = MaoOpto.PosTuning ./ MaxInROIOn; % normalize
MaoOpto.PosTuningNorm(isnan(MaoOpto.PosTuningNorm)) = 0;

%%% CRITERIA #1: place fields must be a continuous region with MinPlaceFieldSize and MinPlaceFieldSize (cm) in width, within which the activity magnitude 
%  must be above 30% of the difference between the maximum and minimum activity in the position tuning curve;

PosTuning = MaoOpto.PosTuning;
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
MaoOpto.PotPlaceFieldPos = PotPlaceFieldPos;
MaoOpto.PotPlaceFieldLength = PotPlaceFieldLength;
MaoOpto.Criteria1Passed = Criteria1Passed;
MaoOpto.Treshold = Treshold;


%%% CRITERIA 2: the mean in-field activity must be at least three times larger than the mean out-of-field activity;  
PosTuning = MaoOpto.PosTuning;
PotPlaceFieldPos = MaoOpto.PotPlaceFieldPos;
PotPlaceFieldLength = MaoOpto.PotPlaceFieldLength;
Criteria1Passed = MaoOpto.Criteria1Passed;

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
MaoOpto.PotPlaceFieldPos2 = MaoOpto.PotPlaceFieldPos;
MaoOpto.PotPlaceFieldPos2(ActRatioInOutPlaceField<InOutRatio) = NaN;
MaoOpto.PotPlaceFieldLength2 = MaoOpto.PotPlaceFieldLength;
MaoOpto.PotPlaceFieldLength2(ActRatioInOutPlaceField<InOutRatio) = NaN;
MaoOpto.Criteria12Passed = Criteria12Passed;


%%% CRITERIA 3: more than one third of the trials must have peak position-mapped activity fall within the potential place field
Reliability = NaN(ROINu,1);
PotPlaceFieldPos2 = MaoOpto.PotPlaceFieldPos2;
PotPlaceFieldLength2 = MaoOpto.PotPlaceFieldLength2;
Criteria12Passed = MaoOpto.Criteria12Passed;
TrialNu = length(Trials);

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

MaoOpto.PlaceFieldStartBin = PotPlaceFieldPos3;
MaoOpto.PlaceFieldBinLength = PotPlaceFieldLength3;
MaoOpto.Criteria123Passed = Criteria123Passed;
MaoOpto.Reliability = Reliability;
MaoOpto.PlaceCells = PlaceCells;

% PLOT place cell ROI averaged (all trial activity)
MaoOpto.placeROINu = length(MaoOpto.PlaceCells);  
MaoOpto.placeROINormActBin = MaoOpto.PosTuningNorm(MaoOpto.PlaceCells,:); % collect only place cell data
makedir = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCellMoreOpto');

if type == 0
    sData.imdata.MaoOpto_dff{Protocol} = struct;
    sData.imdata.MaoOpto_dff{Protocol} = MaoOpto;
    mkdir(makedir,strcat('\Prot',num2str(Protocol),'-dff'));
    SavePath = strcat(makedir,'\Prot',num2str(Protocol),'-dff');
    fn = 'dff';
elseif type == 1
    sData.imdata.MaoOpto_deconv{Protocol} = struct;
    sData.imdata.MaoOpto_deconv{Protocol} = MaoOpto;
    mkdir(makedir,strcat('\Prot',num2str(Protocol),'-deconv'));
    SavePath = strcat(makedir,'\Prot',num2str(Protocol),'-deconv');
    fn = 'deconv';
elseif type == 2
    sData.imdata.MaoOpto_deconv{Protocol} = struct;
    sData.imdata.MaoOpto_deconv{Protocol} = MaoOpto;
    mkdir(makedir,strcat('\Prot',num2str(Protocol),'-spikerate'));
    SavePath = strcat(makedir,'\Prot',num2str(Protocol),'-spikerate');
    fn = 'spikerate';
end
% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% plot position tuning curves with potential place fields
% plot only place cells:
for j = 1:1:length(MaoOpto.PlaceCells)
    i = MaoOpto.PlaceCells(j);
    figure();
    Xaxis = BinSize:BinSize:BinSize*BinNu;
    Ymax = (max(MaoOpto.PosTuning(i,:)))+0.00000000001;
    rectangle('Position',[MaoOpto.PlaceFieldStartBin(i,1)*BinSize Ymax/500 MaoOpto.PlaceFieldBinLength(i,1)*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
    plot(Xaxis,MaoOpto.PosTuning(i,1:BinNu),'LineWidth',2); hold on; % pos tuning curve
    line([0 BinSize*BinNu],[MaoOpto.Treshold(i) MaoOpto.Treshold(i)],'Color','red','LineWidth',1); hold on; 
    line([MaoOpto.PotPlaceFieldPos(i,1)*BinSize MaoOpto.PotPlaceFieldPos(i,1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; line([(MaoOpto.PotPlaceFieldPos(i,1)+MaoOpto.PotPlaceFieldLength(i,1))*BinSize (MaoOpto.PotPlaceFieldPos(i,1)+MaoOpto.PotPlaceFieldLength(i,1))*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
    line([MaoOpto.PotPlaceFieldPos(i,2)*BinSize MaoOpto.PotPlaceFieldPos(i,2)*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on; line([(MaoOpto.PotPlaceFieldPos(i,2)+MaoOpto.PotPlaceFieldLength(i,2))*BinSize (MaoOpto.PotPlaceFieldPos(i,2)+MaoOpto.PotPlaceFieldLength(i,2))*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(i),'-Prot',num2str(Protocol),fn));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(i),'-postun-',fn);
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
end
close all


% sort place cell activity based on peak activity
dataForSortPC = MaoOpto.PosTuning(MaoOpto.PlaceCells,:);
MaxData = max(dataForSortPC, [], 2); % search the max in each row
MaxDataBin = NaN(MaoOpto.placeROINu,2);
MaxDataBin(:,2) = MaoOpto.PlaceCells; % search Max position in bins
for i = 1:1:MaoOpto.placeROINu
    if find(dataForSortPC(i,:)== MaxData(i)) % needed if the max number is twice in the dataset
        MaxDataBin(i,1) = find(dataForSortPC(i,:)== MaxData(i));   % search Max position in bins
        continue
    end
end
SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
MaoOpto.PlaceCellActSortingOrder = SortingOrder(:,2);
%SortedData = NaN(PCNu,BinNu);
SortedData = MaoOpto.PosTuning(SortingOrder(:,2),:); % new matrix containing sorted data.
% Normalize
MaxData2 = max(SortedData, [], 2);
SortedDataMax = SortedData ./ MaxData2;

%PLOT FIGURE NORM
figure('Color','white'); % smoothdata(SortedDataMax,2,'gaussian',5)
imagesc(1:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins,1:MaoOpto.placeROINu,SortedDataMax) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(jet);
c.Label.String = fn;
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 MaoOpto.placeROINu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 MaoOpto.placeROINu],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 MaoOpto.placeROINu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 MaoOpto.placeROINu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 MaoOpto.placeROINu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 MaoOpto.placeROINu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 MaoOpto.placeROINu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 MaoOpto.placeROINu],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat(sData.sessionInfo.fileID,' Mean Place-cell activity - Prot',num2str(Protocol))); 
savefig(fullfile(SavePath,strcat('_',sData.sessionInfo.fileID,'SortedMeanPlaceROI-',fn,'Prot',num2str(Protocol),'.fig')));
saveas(gcf,(fullfile(SavePath,strcat('_',sData.sessionInfo.fileID,'SortedMeanPlaceROI-',fn,'.jpg'))));

end
