function [adata,sData] = learningBehav2(sData,r,adata,savePath)

% BAsed on animals behavior it seemed that they accomodate to a new task quickly usually within 8-10 laps, maximum 15 laps. 
% So I called the behavioral change within the early 1-15 laps accomodation , and later changes development (e.g. between laps 20-30 and 90-100).

% the animal  performs well the task, if:
% 50% of the licks should be above BinStartLickBeforeReward=70 bin. It means, that the animal licks 10 times more/bin in bin 70-80 compared to bin 20-70
% check when it happens in 3-5-10 consecutive trials

% the velocity in the first bin of the RZ should be less than 25 percent of max speed, if min speed in the trial is considered 0
% check when it happens in 3-5-10 consecutive trials

nBins = sData.behavior.meta.nBins;
[nTrials, nBinsPlus] = size(sData.behavior.binning.lickBinned);
RewardZoneSize = 6; % cm
AccomodationLaps = 15;
DevelopmentTestLaps = 10;
VeloAtRZSpeedHigh = 25; % percent, velocity in speed range. Speed range= maxspeed - minspeed in a trial
VeloAtRZSpeedLow = 10;

% generate binned lick data for every trial until reward is given in each trials
waterGivenMatrix = sData.behavior.binning.waterGivenBinned; 
waterGivenMatrix(:,1:10) = 0; % set zero the beginning of trials, I want to see the reward at the end (+10 bin is given at the end)
LickMatrixBeforeRew = sData.behavior.binning.lickBinned;
waterGivenFirst = NaN(nTrials,1);
for i = 1:1:nTrials
    waterGivenFirst(i) = min(find(waterGivenMatrix(i,:) >= 1));
    LickMatrixBeforeRew(i,waterGivenFirst(i):end) = 0; % delete licks after reward given
end
% delete licks before 37 cm, because sometimes animals licks the reward later than RZ, and it would contaminate lick calculation
BinStartLickCalculation = round(37/sData.behavior.meta.binSize); % set 37 because then I will calculate the licks between 37-137 cm = 100 cm
BinStartLickBeforeReward = round(137/sData.behavior.meta.binSize); % expert animal should lick from 20 cm before reward: 137 cm
LickMatrixBeforeRew(:,1:BinStartLickCalculation) = 0;
nBinsEarlyLicks = BinStartLickBeforeReward - 1 - BinStartLickCalculation;
nBinsLateLicks = waterGivenFirst - BinStartLickBeforeReward;

%%% Calculate lick performance
learning2 = struct;
learning2.lick = struct;
learning2.lick.lickEarlyPerBin = sum(LickMatrixBeforeRew(:,BinStartLickCalculation+1:BinStartLickBeforeReward-1),2)/nBinsEarlyLicks;
learning2.lick.lickLatePerBin = NaN(nTrials,1);
for i = 1:1:nTrials
    learning2.lick.lickLatePerBin(i) = sum(LickMatrixBeforeRew(i,BinStartLickBeforeReward:waterGivenFirst(i)))/nBinsLateLicks(i); % Array , collecting each trials lick numbers before reward (137-157 cm) 
end

learning2.lick.LateEarlyRatio = learning2.lick.lickLatePerBin ./ learning2.lick.lickEarlyPerBin;

learning2.lick.details.lickEarlyMean = mean(learning2.lick.lickEarlyPerBin);
learning2.lick.details.lickLateMean = mean(learning2.lick.lickLatePerBin);
learning2.lick.details.lickEarlyStd = std(learning2.lick.lickEarlyPerBin);
learning2.lick.details.lickLateStd = std(learning2.lick.lickLatePerBin);


% position of the 25, 50, 75% quantiles of the lick distibution
learning2.lick.quantilesBins(nTrials,3) = NaN;
learning2.lick.quantilesCm(nTrials,3) = NaN;
for i = 1:1:nTrials
    temp = LickMatrixBeforeRew(i,BinStartLickCalculation:nBins);
    LicksPositionArray = NaN(1,sum(temp));
    counter = 1;
    for j=1:1:nBins-BinStartLickCalculation+1
        if temp(j)>0
            for k = 1:1:temp(j)
                LicksPositionArray(counter) = j + BinStartLickCalculation;
                counter = counter + 1;
            end
        end
    end
    learning2.lick.quantilesBins(i,:) = round(quantile(LicksPositionArray,[0.25 0.50 0.75])); %/nBins*100
    learning2.lick.quantilesCm(i,:) = quantile(LicksPositionArray,[0.25 0.50 0.75])/nBins*sData.behavior.stats.LapLengthCm; %/nBins*100
end
% 50% of the licks should be above BinStartLickBeforeReward=70 bin. It means, that the animal licks 10 times more/bin in bin 70-80 compared to bin 20-70

% In which trial is the 3/5/10th when the animal 50 quantile is at >= as bin 70?
learning2.lick.Q50.TrialsLargerThan70Bin = find(learning2.lick.quantilesBins(:,2) >= BinStartLickBeforeReward);
L = length(learning2.lick.Q50.TrialsLargerThan70Bin);
learning2.lick.Q50.TrialsLargerThan70Bin(L+1:nTrials) = NaN;
learning2.lick.Q50.NotConsTrials3L70 = learning2.lick.Q50.TrialsLargerThan70Bin(3);
learning2.lick.Q50.NotConsTrials5L70 = learning2.lick.Q50.TrialsLargerThan70Bin(5);
learning2.lick.Q50.NotConsTrials10L70 = learning2.lick.Q50.TrialsLargerThan70Bin(10);
% Find the trial when it is true that in 3/5/7 consecutive trials 50Q is larger= than 70bin
counterArray = NaN(10,1);
counter = 0;
learning2.lick.Q50.ConsTrials3L70 = NaN;
learning2.lick.Q50.ConsTrials5L70 = NaN;
learning2.lick.Q50.ConsTrials10L70 = NaN;
for i = 1:1:nTrials
    if learning2.lick.quantilesBins(i,2) >= BinStartLickBeforeReward
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(learning2.lick.Q50.ConsTrials3L70) 
           learning2.lick.Q50.ConsTrials3L70 = i;
       elseif counter == 5 && isnan(learning2.lick.Q50.ConsTrials5L70)
           learning2.lick.Q50.ConsTrials5L70 = i;
       elseif counter == 10 && isnan(learning2.lick.Q50.ConsTrials10L70)
           learning2.lick.Q50.Constrials10L70 = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end
% percentage of trials which is Q50 above 70 bin
learning2.lick.Q50.PercentageOfTrialsL70 = sum(~isnan(learning2.lick.Q50.TrialsLargerThan70Bin))/nTrials;
% as first 8-15 laps usualy needed for accomodation to the new task , calculate percentage of Q50 L70 trials from trial 16
TrialsLarger70BinAfterLap15 = learning2.lick.Q50.TrialsLargerThan70Bin(learning2.lick.Q50.TrialsLargerThan70Bin>AccomodationLaps);
learning2.lick.Q50.PercentageOfTrialsL70AfterT15 = length(TrialsLarger70BinAfterLap15)/(nTrials-AccomodationLaps);
% claculate the change in lap 1-15 : Accomodation (3-3 laps)
learning2.lick.Q50.Accomodation_1_3 = nanmean(learning2.lick.quantilesBins(1:3,2));
learning2.lick.Q50.Accomodation_12_15 = nanmean(learning2.lick.quantilesBins(AccomodationLaps-3:AccomodationLaps,2));
learning2.lick.Q50.AccomodationSbtrBins = learning2.lick.Q50.Accomodation_12_15 - learning2.lick.Q50.Accomodation_1_3;
learning2.lick.Q50.AccomodationRatio = learning2.lick.Q50.Accomodation_12_15 / learning2.lick.Q50.Accomodation_1_3;
% claculate the change after lap 15 : Development (10-10 laps) DevelopmentTestLaps
learning2.lick.Q50.Development_First15 = nanmean(learning2.lick.quantilesBins(AccomodationLaps+1:AccomodationLaps+DevelopmentTestLaps,2));
learning2.lick.Q50.Development_Last15 = nanmean(learning2.lick.quantilesBins(nTrials-DevelopmentTestLaps:nTrials,2));
learning2.lick.Q50.DevelopmentSbtrBins = learning2.lick.Q50.Development_Last15 - learning2.lick.Q50.Development_First15;
learning2.lick.Q50.DevelopmentRatio = learning2.lick.Q50.Development_Last15 / learning2.lick.Q50.Development_First15;


%%% Velocitiy performance
learning2.velo = struct;

[learning2.velo.veloMax, learning2.velo.veloMaxPosBin]= max(sData.behavior.binning.veloBinned(1:nTrials,1:nBinsPlus),[],2);
[learning2.velo.veloMin, learning2.velo.veloMinPosBin] = min(sData.behavior.binning.veloBinned(1:nTrials,1:nBinsPlus),[],2);

% velo at RZ start position (bin1) on max-min scale (0 = the min velo is at RZ, 1= the max velo is there)
learning2.velo.RZspeed.VeloAtRZ = NaN(nTrials,3);
learning2.velo.RZspeed.VeloAtRZ(1:nTrials,1) = (sData.behavior.binning.veloBinned(:,nBins+1) - learning2.velo.veloMin) ./ (learning2.velo.veloMax-learning2.velo.veloMin)*100;
learning2.velo.RZspeed.VeloAtRZ(1:nTrials,2) = (sData.behavior.binning.veloBinned(:,nBins+2) - learning2.velo.veloMin) ./ (learning2.velo.veloMax-learning2.velo.veloMin)*100;
learning2.velo.RZspeed.VeloAtRZ(1:nTrials,3) = (sData.behavior.binning.veloBinned(:,nBins+3) - learning2.velo.veloMin) ./ (learning2.velo.veloMax-learning2.velo.veloMin)*100;

learning2.velo.RZspeed.veloLess25p = find(learning2.velo.RZspeed.VeloAtRZ(:,1) < VeloAtRZSpeedHigh);
L = length(learning2.velo.RZspeed.veloLess25p);
learning2.velo.RZspeed.veloLess25p(L+1:nTrials) = NaN;
learning2.velo.RZspeed.veloLess10p = find(learning2.velo.RZspeed.VeloAtRZ(:,1) < VeloAtRZSpeedLow);
L = length(learning2.velo.RZspeed.veloLess10p);
learning2.velo.RZspeed.veloLess10p(L+1:nTrials) = NaN;

% In which trial is the 3/5/10th when at the RZ the animal speed is below 10/25 percent of max speed?
learning2.velo.RZspeed.NotConsTrials3S25 = learning2.velo.RZspeed.veloLess25p(3);
learning2.velo.RZspeed.NotConsTrials5S25 = learning2.velo.RZspeed.veloLess25p(5);
learning2.velo.RZspeed.NotConsTrials10S25 = learning2.velo.RZspeed.veloLess25p(10);

learning2.velo.RZspeed.NotConsTrials3S10 = learning2.velo.RZspeed.veloLess10p(3);
learning2.velo.RZspeed.NotConsTrials5S10 = learning2.velo.RZspeed.veloLess10p(5);
learning2.velo.RZspeed.NotConsTrials10S10 = learning2.velo.RZspeed.veloLess10p(10);


% Find the trial when it is true that in 3/5/7 consecutive trials the speed is below 10/20 percent of max-min at RZ
% if speed is below 25% at RZ
counterArray = NaN(10,1);
counter = 0;
learning2.velo.RZspeed.ConsTrials3S25 = NaN;
learning2.velo.RZspeed.ConsTrials5S25 = NaN;
learning2.velo.RZspeed.ConsTrials10S25 = NaN;
for i = 1:1:10%nTrials
    if learning2.velo.RZspeed.VeloAtRZ(i,1) < VeloAtRZSpeedHigh
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(learning2.velo.RZspeed.ConsTrials3S25)
           learning2.velo.RZspeed.ConsTrials3S25 = i;
       elseif counter == 5 && isnan(learning2.velo.RZspeed.ConsTrials5S25)
           learning2.velo.RZspeed.ConsTrials5S25 = i;
       elseif counter == 10 && isnan(learning2.velo.RZspeed.ConsTrials10S25)
           learning2.velo.RZspeed.ConsTrials10S25 = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end
% if speed is below 10% at RZ
counterArray = NaN(10,1);
counter = 0;
learning2.velo.RZspeed.ConsTrials3S10 = NaN;
learning2.velo.RZspeed.ConsTrials5S10 = NaN;
learning2.velo.RZspeed.ConsTrials10S10 = NaN;
for i = 1:1:nTrials
    if learning2.velo.RZspeed.VeloAtRZ(i,1) < VeloAtRZSpeedLow
       counter = counter + 1;
       counterArray(counter,1) = i;
       if counter == 3 && isnan(learning2.velo.RZspeed.ConsTrials3S10)
           learning2.velo.RZspeed.ConsTrials3S10 = i;
       elseif counter == 5 && isnan(learning2.velo.RZspeed.ConsTrials5S10)
           learning2.velo.RZspeed.ConsTrials5S10 = i;
       elseif counter == 10 && isnan(learning2.velo.RZspeed.ConsTrials10S10)
           learning2.velo.RZspeed.ConsTrials10S10 = i;
       end
    else 
       counterArray = NaN(10,1);
       counter = 0; 
    end
end


% How many percantage of the trials in this session is the speed below 10/25 percantage of max-min speed?
learning2.velo.RZspeed.TrialsRZSpeedLess25Percent = sum(~isnan(learning2.velo.RZspeed.veloLess25p))/nTrials*100;
learning2.velo.RZspeed.TrialsRZSpeedLess10Percent = sum(~isnan(learning2.velo.RZspeed.veloLess10p))/nTrials*100;

TrialsSpeedLess25AfterLap15 = learning2.velo.RZspeed.veloLess25p(learning2.velo.RZspeed.veloLess25p>AccomodationLaps);
TrialsSpeedLess10AfterLap15 = learning2.velo.RZspeed.veloLess10p(learning2.velo.RZspeed.veloLess10p>AccomodationLaps);
learning2.velo.RZspeed.TrialsRZSpeedLess25PercentAfterT15 = length(TrialsSpeedLess25AfterLap15)/(nTrials-AccomodationLaps)*100;
learning2.velo.RZspeed.TrialsRZSpeedLess10PercentAfterT15 = length(TrialsSpeedLess10AfterLap15)/(nTrials-AccomodationLaps)*100;

% claculate the change in lap 1-15 : Accomodation (3-3 laps)
learning2.velo.RZspeed.Accomodation_1_3 = nanmean(learning2.velo.RZspeed.VeloAtRZ(1:3,1));
learning2.velo.RZspeed.Accomodation_12_15 = nanmean(learning2.velo.RZspeed.VeloAtRZ(AccomodationLaps-3:AccomodationLaps,1));
learning2.velo.RZspeed.AccomodationSbtr = nanmean(learning2.velo.RZspeed.VeloAtRZ(1:3,1)) - nanmean(learning2.velo.RZspeed.VeloAtRZ(AccomodationLaps-3:AccomodationLaps,1));
learning2.velo.RZspeed.AccomodationRatio = nanmean(learning2.velo.RZspeed.VeloAtRZ(1:3,1)) ./ nanmean(learning2.velo.RZspeed.VeloAtRZ(AccomodationLaps-3:AccomodationLaps,1));

% claculate the change after lap 15 : Development (10-10 laps) DevelopmentTestLaps
learning2.velo.RZspeed.Development_First15 = nanmean(learning2.velo.RZspeed.VeloAtRZ(AccomodationLaps+1:AccomodationLaps+DevelopmentTestLaps,1));
learning2.velo.RZspeed.Development_Last15 = nanmean(learning2.velo.RZspeed.VeloAtRZ(nTrials-DevelopmentTestLaps:nTrials,1));
learning2.velo.RZspeed.DevelopmentSbtrBins = learning2.velo.RZspeed.Development_Last15 - learning2.velo.RZspeed.Development_First15;
learning2.velo.RZspeed.DevelopmentRatio = learning2.velo.RZspeed.Development_Last15 / learning2.velo.RZspeed.Development_First15;


% Save file to same path where LV files can be found 
sData.behavior.learning = learning2;
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% for copying data to excel
%r=1;
%adata = NaN(20,36);
list = strsplit(sData.sessionInfo.sessionID,'-');
sessionID = strcat(list(2),list(3));
session = str2double(sessionID)-2019000000;

adata(r,1) = session;
adata(r,2) = sData.behavior.learning.lick.Q50.NotConsTrials3L70;
adata(r,3) = sData.behavior.learning.lick.Q50.NotConsTrials5L70;
adata(r,4) = sData.behavior.learning.lick.Q50.NotConsTrials10L70;
adata(r,5) = sData.behavior.learning.lick.Q50.ConsTrials3L70;
adata(r,6) = sData.behavior.learning.lick.Q50.ConsTrials5L70;
adata(r,7) = sData.behavior.learning.lick.Q50.ConsTrials10L70;
adata(r,8) = sData.behavior.learning.lick.Q50.PercentageOfTrialsL70;
adata(r,9) = sData.behavior.learning.lick.Q50.PercentageOfTrialsL70AfterT15;
adata(r,10) = sData.behavior.learning.lick.Q50.Accomodation_1_3;
adata(r,11) = sData.behavior.learning.lick.Q50.Accomodation_12_15;
adata(r,12) = sData.behavior.learning.lick.Q50.AccomodationSbtrBins;
adata(r,13) = sData.behavior.learning.lick.Q50.AccomodationRatio;
adata(r,14) = sData.behavior.learning.lick.Q50.Development_First15;
adata(r,15) = sData.behavior.learning.lick.Q50.Development_Last15;
adata(r,16) = sData.behavior.learning.lick.Q50.DevelopmentSbtrBins;
adata(r,17) = sData.behavior.learning.lick.Q50.DevelopmentRatio;
adata(r,18) = NaN;
adata(r,19) = sData.behavior.learning.velo.RZspeed.NotConsTrials3S25;
adata(r,20) = sData.behavior.learning.velo.RZspeed.NotConsTrials5S25;
adata(r,21) = sData.behavior.learning.velo.RZspeed.NotConsTrials10S25;
adata(r,22) = sData.behavior.learning.velo.RZspeed.NotConsTrials3S10;  
adata(r,23) = sData.behavior.learning.velo.RZspeed.NotConsTrials5S10; 
adata(r,24) = sData.behavior.learning.velo.RZspeed.NotConsTrials10S10; 
adata(r,25) = sData.behavior.learning.velo.RZspeed.ConsTrials3S25;
adata(r,26) = sData.behavior.learning.velo.RZspeed.ConsTrials5S25;
adata(r,27) = sData.behavior.learning.velo.RZspeed.ConsTrials10S25;
adata(r,28) = sData.behavior.learning.velo.RZspeed.TrialsRZSpeedLess25Percent;
adata(r,29) = sData.behavior.learning.velo.RZspeed.TrialsRZSpeedLess10Percent;
adata(r,30) = sData.behavior.learning.velo.RZspeed.TrialsRZSpeedLess25PercentAfterT15;
adata(r,31) = sData.behavior.learning.velo.RZspeed.TrialsRZSpeedLess10PercentAfterT15;
adata(r,32) = sData.behavior.learning.velo.RZspeed.Accomodation_1_3; 
adata(r,33) = sData.behavior.learning.velo.RZspeed.Accomodation_12_15;
adata(r,34) = sData.behavior.learning.velo.RZspeed.AccomodationSbtr;
adata(r,35) = sData.behavior.learning.velo.RZspeed.Development_First15;
adata(r,36) = sData.behavior.learning.velo.RZspeed.Development_Last15;
adata(r,37) = sData.behavior.learning.velo.RZspeed.DevelopmentSbtrBins;
adata(r,38) = sData.behavior.learning.velo.RZspeed.DevelopmentRatio;

end