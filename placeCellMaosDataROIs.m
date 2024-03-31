function sData = placeCellMaosDataROIs(sData,type,InOutRatio,ReliabiliyIndex,ROIsDiscard)
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s (set in calcCADATAVelMinNori function input). 
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)
% My data is already filtered for low activity in CaPreProc function(the 99% percentile of ROI data should have been > 3*SD of noise of ROI data)

%%%PARAMATERS TO SET:
%type = 0; %0 use dff data, 1 use deconvolved data, 2 use firing rate data, -1: use dff slow removed data
% InOutRatio = 3; Reliability = 0.34; regarding Mao, for VIP cells IOR = 2, Rel=o.2
%IsOpto = 0; % 0 if it is not and optical session, 1 if it is 
%ROIsDiscard: array containing ROI IDs to discard. 

MinPlaceFieldSize = 4; % cm minimum size of place field Mao:15
MaxPlaceFieldSize = 80; % cm maximum size of place field Mao:120
%CRITERIA 1 the maximum activity in the potential place field must have been at least 1+activityTreshold -fold larger (Mao used 1.3x for deconv signal) compared to minimum activity in pos tuning curve
if type == 0 || type == -1 % for dff data
    activityTreshold = 0.6;  
else % deconvolved signal
    activityTreshold = 0.3; % Mao = 0.3
    GaussianFilter = round(sData.behavior.meta.binSize*5); % I set 5 bins to smooth (10 cm) . 2 cm Binsize usually
end
%CRITERIA 2 the mean in-field activity must be at least three times larger than the mean out-of-field activity;
%InOutRatio = 2; % Mao = 3
%CRITERIA 3 more than one third of the trials must have peak position-mapped activity fall within the potential place field
%ReliabiliyIndex = 0.2; % Mao = 0.34

%%% I set place cell tuning around peak to be 60% as larger as minimum, and place field act should be 2x larger than out of field
%%% INITIALIZE PARAMETERS 
% Load parameters:
Mao = struct; % collect place cell data
Mao.params.activityTreshold = activityTreshold;
Mao.params.InOutRatio = InOutRatio;
Mao.params.ReliabiliyIndex = ReliabiliyIndex;
BinNu = sData.behavior.meta.nBins;
TRNu = sData.behavior.wheelLapImaging;
BinSize = sData.behavior.meta.binSize;
ROINu = sData.imdata.nROIs;


if type == 0 
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-dFF-selection');
    SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-dFF-selection');
elseif type == -1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-dFF-SR-selection');
    SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-dFF-SR-selection');
elseif type == 1
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-deconv-selection');
    SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-deconv-selection');
elseif type == 2
    mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'MaoPlaceCells-FR-selection');
    SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\MaoPlaceCells-FR-selection');
end

% Create arrays for processed data:
dataset = struct;
if type == 0
    dataset = sData.imdata.binned.RoidFF;
    Mao.datatype = 'dFF';
elseif type == -1
    dataset = sData.imdata.binned.RoidFFSR;
    Mao.datatype = 'dFF-slowremoved';
elseif type == 1
    dataset = sData.imdata.binned.RoiDeconvolved;
    Mao.datatype = 'deconvolved';
elseif type == 2
    dataset = sData.imdata.binned.RoiSpikeRate;
    Mao.datatype = 'spikerate';
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

Mao.PosTuning = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    for j = 1:1:BinNu
        Mao.PosTuning(i,j) = nanmean(dataset{1,i}(:,j)); %dff or deocnvolved or SpikeRate        
    end
end

%%% USE 3-BIN GAUSSIAN SMOOOTHING for deconvolved data.
if type == 1
    GaussianSmoothed = NaN(TRNu,BinNu);
    for i = 1:1:ROINu
        GaussianSmoothed(i,:) = smoothdata(Mao.PosTuning(i,:),'gaussian',GaussianFilter);
    end
    Mao.PosTuning = GaussianSmoothed;
end


% Normalize mean activity for visualization
Mao.PosTuningNorm = NaN(ROINu,BinNu);
MaxInROI = max(Mao.PosTuning,[],2); % search for maximum value (in a bin) in each ROI
Mao.PosTuningNorm = Mao.PosTuning ./ MaxInROI; % normalize

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
MinInROI = min(Mao.PosTuning,[],2);
MinInROI(MinInROI<0) = 0;
Treshold = MinInROI + (MaxInROI - MinInROI)*activityTreshold;   %%% CHANGE 0.3 to 0.6  % activity must be larger in place field than this value (130% of minimum)
LargerThanTresholdBins = zeros(ROINu,BinNu);
for i = 1:1:ROINu
    SignalTemp3 = Mao.PosTuning(i,:);
    SignalTemp3(SignalTemp3<Treshold(i))=NaN; % put NaN to bins where activity is smaller than treshold
    LargerThanTresholdBins(i,:) = SignalTemp3; % matrix tresholded activity        
end
LargerThanTresholdBins2 = ~isnan(LargerThanTresholdBins); % ones when the activity is above treshold. Count ones next to each other
%SignalTemp4 = circshift(LargerThan30Bins2,-1,1); % copy matrix, and shift up to lenghten trial to check place field size if next trial is needed
%SignalTemp4(ROINu,:) = zeros;
ExtraZeros = zeros(ROINu,1);
LargerThanTresholdBins3 = [ExtraZeros LargerThanTresholdBins2 LargerThanTresholdBins2]; % concatenate three matrices to lengten trials 
LargerThanTresholdBins4 = NaN(ROINu,2*BinNu+1); % calculate potential place field size. Min: 4 cm, max: 100 cm.
kmax = MaxPlaceFieldSize/BinSize; % max place field is 100 cm, do not check biger size
for i = 1:1:ROINu
    for j = 1:1:2*BinNu % + kmax is because next trial check
      if LargerThanTresholdBins3(i,j+1) == 1 % find the first / next one, in this matrix data are shifted with one column ->(j+1)
          for k = 1:1:BinNu+kmax-j % count how many ones are after it
              if LargerThanTresholdBins3(i,j+1+k) == 0 % if there is a zero (low activity), write the length of place field and jump to the next ones           
                  LargerThanTresholdBins4(i,j+1) = k; % how many bin long place field is it , I put an extra column in the beginning
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
Mao.Criteria1Passed = NaN(ROINu,1);
FieldMinBin = round(MinPlaceFieldSize/BinSize);
FieldMaxBin = round(MaxPlaceFieldSize/BinSize);
for i = 1:1:ROINu
    SignalTemp6 = PlaceFieldSize(i,:);
    if any(SignalTemp6<FieldMaxBin) && any(SignalTemp6>FieldMinBin) % at binsize=2 
        Mao.Criteria1Passed(i)=i;
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
Mao.PotPlaceFieldPos = PotPlaceFieldPos;
Mao.PotPlaceFieldLength = PotPlaceFieldLength;


%%% CRITERIA 2: the mean in-field activity must be at least three times larger than the mean out-of-field activity; 
ActAllTrials = [Mao.PosTuning Mao.PosTuning]; % concatenate two matrices to have to trials consecutive
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
Mao.Criteria12Passed = Mao.Criteria1Passed;
for i = 1:1:ROINu
    if ~any(ActRatioInOutPlaceField(i,:) > InOutRatio) %%% CHANGE to 2 - 3?
       Mao.Criteria12Passed(i) = NaN;
    end
end
% update potential place field matrix 
PotPlaceFieldPos2 = PotPlaceFieldPos;
PotPlaceFieldPos2(ActRatioInOutPlaceField<InOutRatio) = NaN;
PotPlaceFieldLength2 = PotPlaceFieldLength;
PotPlaceFieldLength2(ActRatioInOutPlaceField<InOutRatio) = NaN;

%%% CRITERIA 3: more than one third of the trials must have peak position-mapped activity fall within the potential place field
Reliability = NaN(ROINu,1);
Column = NaN(ROINu,1);
for i = 1:1:ROINu
    if any(~isnan(PotPlaceFieldPos2(i,:))) % if no potential place field, jump to next roi
        counter = 0; % counter for reliable trials
        Column(i) = min(find(~isnan(PotPlaceFieldPos2(i,:)))); % if two peaks pass the criteria, choose the first
        PFStart = PotPlaceFieldPos2(i,Column(i)); % place field start (bin)
        PFEnd = PotPlaceFieldPos2(i,Column(i)) + PotPlaceFieldLength2(i,Column(i)); % place field end
        for j = 1:1:TRNu
            TrialData = dataset{1,i}(j,:); % put each trial data into a container
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
        Reliability(i,1) = counter / TRNu;  % reliability data, 1 is the max. Needed >0.33 to be a place cell
    end
end
Mao.Criteria123Passed = Mao.Criteria12Passed;
for i = 1:1:ROINu
    if Reliability(i,1) < ReliabiliyIndex
        Mao.Criteria123Passed(i) = NaN;
    end
end
Mao.Criteria123Passed2 = Mao.Criteria123Passed;
Mao.Criteria123Passed2(ROIsDiscard) = NaN;
Mao.PlaceCells =  Mao.Criteria123Passed2(~isnan(Mao.Criteria123Passed2));
% update potential place field matrix 
PotPlaceFieldPos3 = PotPlaceFieldPos2;
PotPlaceFieldPos3(Reliability < ReliabiliyIndex) = NaN;
PotPlaceFieldLength3 = PotPlaceFieldLength2;
PotPlaceFieldLength3(Reliability < ReliabiliyIndex) = NaN;
Mao.PlaceFieldStartBin = PotPlaceFieldPos3;
Mao.PlaceFieldBinLength = PotPlaceFieldLength3;
Mao.Reliability= Reliability;

% PLOT place cell ROI averaged (all trial activity)
Mao.placeROINu = length(Mao.PlaceCells);  
Mao.placeROINormActBin = Mao.PosTuningNorm(Mao.PlaceCells,:); % collect only place cell data


if type == 0 
    sData.imdata.MaoPC_dff = struct;
    sData.imdata.MaoPC_dff = Mao;
    fn = 'dff';
    dataForSortPC = sData.imdata.binned.MeanRoiAct(Mao.PlaceCells,:);
    dataForSort = sData.imdata.binned.MeanRoiAct;
elseif type == -1 
    sData.imdata.MaoPC_dffSR = struct;
    sData.imdata.MaoPC_dffSR = Mao;
    fn = 'dffSR';
    dataForSortPC = sData.imdata.binned.MeanRoiActSR(Mao.PlaceCells,:);
    dataForSort = sData.imdata.binned.MeanRoiActSR;
elseif type == 1 
    sData.imdata.MaoPC_deconv = struct;
    sData.imdata.MaoPC_deconv = Mao;
    fn = 'deconv';
    dataForSortPC = sData.imdata.binned.MeanRoiAct_Deconv(Mao.PlaceCells,:);
    dataForSort = sData.imdata.binned.MeanRoiAct_Deconv;
elseif type == 2 
    sData.imdata.MaoPC_FR = struct;
    sData.imdata.MaoPC_FR = Mao;
    fn = 'firing rate';
    dataForSortPC = sData.imdata.binned.MeanRoiAct_FR(Mao.PlaceCells,:);
    dataForSort = sData.imdata.binned.MeanRoiAct_FR;
end


% plot position tuning curves with potential place fields
% plot only place cells (in light off trials):
for j = 1:1:length(Mao.PlaceCells)
    i = Mao.PlaceCells(j);
    figure();
    Xaxis = BinSize:BinSize:BinSize*BinNu;
    Ymax = (max(Mao.PosTuning(i,:)))+0.00000000001;
    rectangle('Position',[Mao.PlaceFieldStartBin(i,Column(i))*BinSize Ymax/500 Mao.PlaceFieldBinLength(i,Column(i))*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
    plot(Xaxis,Mao.PosTuning(i,1:BinNu),'LineWidth',2); hold on; % pos tuning curve
    line([0 BinSize*BinNu],[Treshold(i) Treshold(i)],'Color','red','LineWidth',1); hold on; 
    line([Mao.PotPlaceFieldPos(i,1)*BinSize Mao.PotPlaceFieldPos(i,1)*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on; line([(Mao.PotPlaceFieldPos(i,1)+Mao.PotPlaceFieldLength(i,1))*BinSize (Mao.PotPlaceFieldPos(i,1)+Mao.PotPlaceFieldLength(i,1))*BinSize],[0 Ymax],'Color',[1 0.3 0.2],'LineStyle','--','LineWidth',1); hold on;
    line([Mao.PotPlaceFieldPos(i,2)*BinSize Mao.PotPlaceFieldPos(i,2)*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on; line([(Mao.PotPlaceFieldPos(i,2)+Mao.PotPlaceFieldLength(i,2))*BinSize (Mao.PotPlaceFieldPos(i,2)+Mao.PotPlaceFieldLength(i,2))*BinSize],[0 Ymax],'Color',[0.6 0.3 0],'LineStyle','--','LineWidth',1); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(i),'-',fn));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(i),fn);
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
 end
close all


% sort place cell activity based on peak activity
PCNu = length(Mao.PlaceCells);
MaxData = max(dataForSortPC, [], 2); % search the max in each row
MaxDataBin = NaN(PCNu,2);
MaxDataBin(:,2) = Mao.PlaceCells; % search Max position in bins
for i = 1:1:PCNu
    if find(dataForSortPC(i,:)== MaxData(i)) % needed if the max number is twice in the dataset
        MaxDataBin(i,1) = find(dataForSortPC(i,:)== MaxData(i));   % search Max position in bins
        continue
    end
end
SortingOrder = sortrows(MaxDataBin, 1); % second column is the sorted ROI order. plot in these order
Mao.PlaceCellActSortingOrder = SortingOrder(:,2);
%SortedData = NaN(PCNu,BinNu);
SortedData = dataForSort(SortingOrder(:,2),:); % new matrix containing sorted data.
% Normalize
MaxData2 = max(SortedData, [], 2);
SortedDataMax = SortedData ./ MaxData2;

%PLOT FIGURE NORM
figure('Color','white'); % smoothdata(SortedDataMax,2,'gaussian',5)
imagesc(1:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins,1:PCNu,smoothdata(SortedDataMax,2,'gaussian',5)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(jet);
c.Label.String = strcat(fn,'(normalized)');
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
%ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
%set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(ax,'yticklabel',[])
title(strcat(sData.sessionInfo.fileID,'Place-cell activity sorted based on peak activity')); 
if ~isempty(Mao.PlaceCellActSortingOrder)
    text(-20,(length(Mao.PlaceCellActSortingOrder)/2)+0.5,num2str(Mao.PlaceCellActSortingOrder),'FontSize',round(200/length(Mao.PlaceCellActSortingOrder)));
end
savefig(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROI-',fn,'.fig')));
saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROI-',fn,'.jpg'))));

%PLOT FIGURE NOT NORM
figure('Color','white'); 
imagesc(1:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins,1:PCNu,smoothdata(SortedData,2,'gaussian',5)) %(1:number of bins;1:number of trials)
c = colorbar;
colormap(jet);
c.Label.String = fn;
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
if type == 0
    caxis([0 0.5]);
else
    caxis([0 1]);
end
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
%ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
%set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(ax,'yticklabel',[])
title(strcat(sData.sessionInfo.fileID,'Place-cell activity sorted based on peak activity')); 
if ~isempty(Mao.PlaceCellActSortingOrder)
    text(-20,(length(Mao.PlaceCellActSortingOrder)/2)+0.5,num2str(Mao.PlaceCellActSortingOrder),'FontSize',round(200/length(Mao.PlaceCellActSortingOrder)));
end
savefig(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROINotNorm.fig')));
saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROINotNorm.jpg'))));

%PLOT FIGURE NOT NORM ZSCORED
figure('Color','white'); 
imagesc(1:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins,1:PCNu,zscore(smoothdata(SortedData,2,'gaussian',5)')') %(1:number of bins;1:number of trials)
c = colorbar;
colormap(jet);
c.Label.String = strcat(fn,'(ZScored)');
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C1B sData.imdata.cues.C1B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C2B sData.imdata.cues.C2B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C3B sData.imdata.cues.C3B],[0 PCNu],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on; line([sData.imdata.cues.C4B sData.imdata.cues.C4B],[0 PCNu],'Color','black','LineStyle','--','LineWidth',2); hold on;
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
Xticks = linspace(1, 160, numel(xticklabels)); set(gca, 'XTick', Xticks, 'XTickLabel', xticklabels)
%ylabel('ROIs'); yticklabels = 0:10:ROINu; yticks = linspace(1, ROINu, numel(yticklabels));
%set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
set(ax,'yticklabel',[])
title(strcat(sData.sessionInfo.fileID,' Place-cell activity sorted based on peak activity')); 
if ~isempty(Mao.PlaceCellActSortingOrder)
    text(-20,(length(Mao.PlaceCellActSortingOrder)/2)+0.5,num2str(Mao.PlaceCellActSortingOrder),'FontSize',round(200/length(Mao.PlaceCellActSortingOrder)));
end
savefig(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROINotNormZScore.fig')));
saveas(gcf,(fullfile(SavePath,strcat(sData.sessionInfo.fileID,'SortedMeanPlaceROINotNormZScore.jpg'))));

% Save file to same path where other files can be found 
if type == 0 
    sData.imdata.MaoPC_dff_selection = Mao; 
elseif type == -1 
    sData.imdata.MaoPC_dffSR_selection = Mao;
elseif type == 1 
    sData.imdata.MaoPC_deconv_selection = Mao;    
elseif type == 2
    sData.imdata.MaoPC_spikerate = Mao;    
end
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end
