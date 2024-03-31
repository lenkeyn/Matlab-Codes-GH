function PC_Mao1 = placecellMao2fordff(sData)
% Calculating which ROI is considered as a place cell based on Mao et al., 2017, criteria#1
% Note: Mao used 1.5 cm bins, I use 2 cm bins and minimum velocity of 1 cm/s (set in calcCADATAVelMinNori function input). 
% Mao summated the activity within a bin and normalized by the occupancy (the time spent) in each bin. I used mean activity within a bin, which is the same (just have to divide all data with SampleTime= 32.3 msec)
% My data is already filtered for low activity in CaPreProc function(the 99% percentile of ROI data should have been > 3*SD of noise of ROI data)


%%% I set place cell tuning around peak to be 60% as larger as minimum, and place field act should be 2x larger than out of field
%%% INITIALIZE PARAMETERS 
% Load parameters:
BinNu = sData.behavior.meta.nBins;
TRNu = sData.behavior.wheelLapImaging;
SampleNu = sData.imdata.nSamples;
FRHz = 31;
VelMin = 0.1;
BinSize = sData.behavior.meta.binSize;
ROINu = sData.imdata.nROIs;
SelectedROIs = (1:1:ROINu); % the row number is the new name of the ROI, the value is the orignal mane in the recording
SavePath = 'C:\MATLAB\SAVE';

% Create arrays for processed data:
PC_Mao1 = struct;
%PC_Mao1.dFF = CADATA.CaBinned.SingleRoidFF; % Mao used deconvolved data
%PC_Mao1.SpikeRate = CADATA.CaBinned.SingleRoiSpikeRate; % SpikeRate is the integrated deconvolved data

%%% USE 3-BIN GAUSSIAN SMOOOTHING. (Mao used 4.5 cm = 3 bin, I will use 3 bin = 6 cm if BinSize is 2 cm)
% Finally not used, since during SpikeRate estimation there is a Gaussian smoothing with FrameRate/2 (15sec), and this new filtering does not cause any effect
%{
for i = 1:1:ROINu
    PC_Mao1.Gauss{i} = NaN(TRNu,BinNu); 
    for j = 1:1:TRNu
        SignalTemp1 = PC_Mao1.SpikeRate{i}(j,:);
        GaussianFiltered = smoothdata(SignalTemp1,'gaussian',3);
        PC_Mao1.Gauss{i}(j,:) = GaussianFiltered;
    end
end
% compare ROI data with and without this Gaussian smoothing:
%roi = 16; % set actual roi 
%figure();
%plotHeatBin(PC_Mao1.Deconv{roi},CADATA.FileID,roi,'spikerate',160,TRNu); hold on;
%plotHeatBin(PC_Mao1.Gauss{roi},CADATA.FilePath,roi,'Gauss-spikerate',160,TRNu); hold on;
%caxis([0 inf]);
%}

%%% GENERATE POSITION TUNING CURVES FOR EACH ROIs
% averaging the position-mapped activity across all trials.
PC_Mao1.dFF_AllTrials = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    for j = 1:1:BinNu
        SignalTemp2 = sData.imdata.binned.RoidFF{1,i}(:,j); %SpikeRate
        PC_Mao1.dFF_AllTrials(i,j) = nanmean(SignalTemp2);
    end
end

% Normalize mean activity for visualization
PC_Mao1.dFF_AllTrialsNorm = NaN(ROINu,BinNu);
MaxInROI = max(PC_Mao1.dFF_AllTrials,[],2); % search for maximum value (in a bin) in each ROI
PC_Mao1.dFF_AllTrialsNorm = PC_Mao1.dFF_AllTrials ./ MaxInROI; % normalize

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

%%% CRITERIA #1: place fields must be a continuous region with 15-120 cm in width, within which the activity magnitude 
%  must be above 30% of the difference between the maximum and minimum activity in the position tuning curve;
MinInROI = min(PC_Mao1.dFF_AllTrials,[],2);
Percent30 = MinInROI + (MaxInROI - MinInROI)*0.6;   %%% CHANGE 0.3 to 0.6  % activity must be larger in place field than this value (130% of minimum)
LargerThan30Bins = zeros(ROINu,BinNu);
for i = 1:1:ROINu
    SignalTemp3 = PC_Mao1.dFF_AllTrials(i,:);
    SignalTemp3(SignalTemp3<Percent30(i))=NaN; % put NaN to bins where activity is smaller than treshold
    LargerThan30Bins(i,:) = SignalTemp3; % matrix tresholded activity        
end
LargerThan30Bins2 = ~isnan(LargerThan30Bins); % ones when the activity is above treshold. Count ones next to each other
%SignalTemp4 = circshift(LargerThan30Bins2,-1,1); % copy matrix, and shift up to lenghten trial to check place field size if next trial is needed
%SignalTemp4(ROINu,:) = zeros;
ExtraZeros = zeros(ROINu,1);
LargerThan30Bins3 = [ExtraZeros LargerThan30Bins2 LargerThan30Bins2]; % concatenate three matrices to lengten trials 
LargerThan30Bins4 = NaN(ROINu,2*BinNu+1); % calculate potential place field size. Min: 16 cm, max: 120 cm.
kmax = 120/BinSize; % max place field is 120 cm, do not check biger size
for i = 1:1:ROINu
    for j = 1:1:2*BinNu % + kmax is because next trial check
      if LargerThan30Bins3(i,j+1) == 1 % find the first / next one, in this matrix data are shifted with one column ->(j+1)
          for k = 1:1:BinNu+kmax-j % count how many ones are after it
              if LargerThan30Bins3(i,j+1+k) == 0 % if there is a zero (low activity), write the length of place field and jump to the next ones           
                  LargerThan30Bins4(i,j+1) = k; % how many bin long place field is it , I put an extra column in the beginning
                  break
              end
          end
      end
   end
end
PlaceFieldSize = NaN(ROINu,BinNu);
for i = 1:1:ROINu
    SignalTemp5 = LargerThan30Bins4(i,1:BinNu+1);
    for j = 1:1:BinNu
        if isnan(SignalTemp5(j)) && ~isnan(SignalTemp5(j+1)) 
            PlaceFieldSize(i,j) = SignalTemp5(j+1); % j-1, because it ws sifted in LargerThan30Bins matrix to the right
        end
    end
end
% Collect ROIs having at least one potential place field
PC_Mao1.Criteria1Passed = NaN(ROINu,1);
FieldMinBin = round(15/BinSize);
FieldMaxBin = round(120/BinSize);
for i = 1:1:ROINu
    SignalTemp6 = PlaceFieldSize(i,:);
    if any(SignalTemp6<FieldMaxBin) && any(SignalTemp6>FieldMinBin) % at binsize=2 cm it is 8 and 60 bin
        PC_Mao1.Criteria1Passed(i)=i;
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
PC_Mao1.PotPlaceFieldPos = PotPlaceFieldPos;
PC_Mao1.PotPlaceFieldLength = PotPlaceFieldLength;


%%% CRITERIA 2: the mean in-field activity must be at least three times larger than the mean out-of-field activity; 
% I check for each potential place field individually (even when there are two potential place fields).
%PotentialPC2 = PC_Mao1.Criteria1Passed(~isnan(PC_Mao1.Criteria1Passed));
ActAllTrials = [PC_Mao1.dFF_AllTrials PC_Mao1.dFF_AllTrials]; % concatenate two matrices to have to trials consecutive
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
PC_Mao1.Criteria12Passed = PC_Mao1.Criteria1Passed;
for i = 1:1:ROINu
    if ~any(ActRatioInOutPlaceField(i,:) > 2) %%% CHANGE 3 to 2
       PC_Mao1.Criteria12Passed(i) = NaN;
    end
end
% update potential place field matrix 
PotPlaceFieldPos2 = PotPlaceFieldPos;
PotPlaceFieldPos2(ActRatioInOutPlaceField<3) = NaN;
PotPlaceFieldLength2 = PotPlaceFieldLength;
PotPlaceFieldLength2(ActRatioInOutPlaceField<3) = NaN;

%%% CRITERIA 3: more than one third of the trials must have peak position-mapped activity fall within the potential place field
Reliability = NaN(ROINu,1);
RI = 0.34; % reliablity index, reliability should be larger in order to be a place cell
for i = 1:1:ROINu
    if any(~isnan(PotPlaceFieldPos2(i,:))) % if no potential place field, jump to next roi
        counter = 0; % counter for reliable trials
        PFStart = PotPlaceFieldPos2(i); % place field start (bin)
        PFEnd = PotPlaceFieldPos2(i) + PotPlaceFieldLength(i); % place field end
        for j = 1:1:TRNu
            TrialData = sData.imdata.binned.RoidFF{1,i}(j,:); % put each trial data into a container
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
PC_Mao1.Criteria123Passed = PC_Mao1.Criteria12Passed;
for i = 1:1:ROINu
    if Reliability(i,1) < RI
        PC_Mao1.Criteria123Passed(i) = NaN;
    end
end
PC_Mao1.PlaceCells =  PC_Mao1.Criteria123Passed(~isnan(PC_Mao1.Criteria123Passed));
% update potential place field matrix 
PotPlaceFieldPos3 = PotPlaceFieldPos2;
PotPlaceFieldPos3(Reliability(i,1) < RI) = NaN;
PotPlaceFieldLength3 = PotPlaceFieldLength2;
PotPlaceFieldLength3(Reliability(i,1) < RI) = NaN;
PC_Mao1.PlaceFieldStartBin = PotPlaceFieldPos3;
PC_Mao1.PlaceFieldBinLength = PotPlaceFieldLength3;

% PLOT place cell ROI averaged (all trial activity)
placeROINu = length(PC_Mao1.PlaceCells);  
PC_Mao1.placeROINormActBin = PC_Mao1.dFF_AllTrialsNorm(PC_Mao1.PlaceCells,:); % collect only place cell data

plotSortROIsMaxActHeat(PC_Mao1.placeROINormActBin);
savefig(fullfile(SavePath,'MeanPlaceROIactNormSorted.fig'));
saveas(gcf,(fullfile(SavePath,'MeanPlaceROIactNormSorted.jpg')));

sData.imdata.MaoPC = PC_Mao1;

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end
