function [copy,sData] = learningBehav(sData)

nBins = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging-1;
LapsAnalysed = 20;
BinsToCollect = round(20/sData.behavior.meta.binSize);
RewardZoneSize = 6; % cm


%%% Lick performance
LickLateStart = nBins-BinsToCollect;
LickMiddleStart = floor(100/sData.behavior.meta.binSize);
LickEarlyStart = floor(60/sData.behavior.meta.binSize);

sData.learning = struct;
sData.learning.lick = struct;
sData.learning.lick.MeanLickSessionPerBin = mean2(sData.behavior.binning.lickBinned);
sData.learning.lick.lickLate = sum(sData.behavior.binning.lickBinned(:,LickLateStart:nBins),2); % Array , collecting each trials lick numbers before reward (137-157 cm) 
sData.learning.lick.lickMiddle = sum(sData.behavior.binning.lickBinned(:,LickMiddleStart:LickMiddleStart+BinsToCollect),2); % Array , collecting each trials lick numbers in the middle of the trial (90-110 cm) 
sData.learning.lick.lickEarly  = sum(sData.behavior.binning.lickBinned(:,LickEarlyStart:LickEarlyStart+BinsToCollect),2);

sData.learning.lick.lickLateMean = mean(sData.learning.lick.lickLate);
sData.learning.lick.lickMiddleMean = mean(sData.learning.lick.lickMiddle);
sData.learning.lick.lickEarlyMean = mean(sData.learning.lick.lickEarly);

sData.learning.lick.lickLateStd = std(sData.learning.lick.lickLate);
sData.learning.lick.lickMiddleStd = std(sData.learning.lick.lickMiddle);
sData.learning.lick.lickEarlyStd = std(sData.learning.lick.lickEarly);

sData.learning.lick.lickDiffLateEarly = sData.learning.lick.lickLate - sData.learning.lick.lickEarly;
sData.learning.lick.lickDiffLateMiddle = sData.learning.lick.lickLate - sData.learning.lick.lickMiddle;
sData.learning.lick.lickDiffMiddleEarly = sData.learning.lick.lickMiddle - sData.learning.lick.lickEarly;

sData.learning.lick.lickDiffLateEarlyMean = mean(sData.learning.lick.lickDiffLateEarly);
sData.learning.lick.lickDiffLateMiddleMean = mean(sData.learning.lick.lickDiffLateMiddle);
sData.learning.lick.lickDiffMiddleEarlyMean = mean(sData.learning.lick.lickDiffMiddleEarly);

sData.learning.lick.lickDiffLateEarlyStd = std(sData.learning.lick.lickDiffLateEarly);
sData.learning.lick.lickDiffLateMiddleStd = std(sData.learning.lick.lickDiffLateMiddle);
sData.learning.lick.lickDiffMiddleEarlyStd = std(sData.learning.lick.lickDiffMiddleEarly);

sData.learning.lick.lickDiffLateEarlyReliability = sData.learning.lick.lickDiffLateEarlyMean / sData.learning.lick.lickDiffLateEarlyStd;
sData.learning.lick.lickDiffLateMiddleMeanReliability = sData.learning.lick.lickDiffLateMiddleMean / sData.learning.lick.lickDiffLateMiddleStd;
sData.learning.lick.lickDiffMiddleEarlyMeanReliability = sData.learning.lick.lickDiffMiddleEarlyMean / sData.learning.lick.lickDiffMiddleEarlyStd;

% Performance
sData.learning.lick.LickPerformance_LateMiddle = sData.learning.lick.lickDiffLateMiddleMean; 
sData.learning.lick.LickPerformance_LateMiddleNorm =  sData.learning.lick.lickDiffLateMiddleMean / sData.learning.lick.lickDiffMiddleEarlyMean;
sData.learning.lick.LickPerformance_LateEarly = sData.learning.lick.lickDiffLateEarlyMean;
sData.learning.lick.LickPerformance_LateEarlyNorm =  sData.learning.lick.lickDiffLateEarlyMean / sData.learning.lick.lickDiffMiddleEarlyMean;


% First and Second half od the recording
sData.learning.lick.FirstSecondHalf = struct;
sData.learning.lick.FirstSecondHalf.lickLateMeanFirst = mean(sData.learning.lick.lickLate(1:LapsAnalysed));
sData.learning.lick.FirstSecondHalf.lickLateMeanLast = mean(sData.learning.lick.lickLate(nTrials-LapsAnalysed:nTrials));
sData.learning.lick.FirstSecondHalf.lickMiddleMeanFirst = mean(sData.learning.lick.lickMiddle(1:LapsAnalysed));
sData.learning.lick.FirstSecondHalf.lickMiddleMeanLast = mean(sData.learning.lick.lickMiddle(nTrials-LapsAnalysed:nTrials));
sData.learning.lick.FirstSecondHalf.lickEarlyMeanFirst = mean(sData.learning.lick.lickEarly(1:LapsAnalysed));
sData.learning.lick.FirstSecondHalf.lickEarlyMeanLast = mean(sData.learning.lick.lickEarly(nTrials-LapsAnalysed:nTrials));
sData.learning.lick.FirstSecondHalf.Development = mean([(sData.learning.lick.FirstSecondHalf.lickLateMeanLast - sData.learning.lick.FirstSecondHalf.lickLateMeanFirst) (sData.learning.lick.FirstSecondHalf.lickMiddleMeanFirst-sData.learning.lick.FirstSecondHalf.lickMiddleMeanLast) (sData.learning.lick.FirstSecondHalf.lickEarlyMeanFirst-sData.learning.lick.FirstSecondHalf.lickEarlyMeanLast)]);


%%% Velocitiy performance
sData.learning.velo = struct;
sData.learning.velo.veloMeanMax = max(sData.behavior.binning.meanVeloBinned(1:nBins));
sData.learning.velo.veloMeanMin = min(sData.behavior.binning.meanVeloBinned(1:nBins));
sData.learning.velo.veloMeanMaxPosBin = find(sData.behavior.binning.meanVeloBinned==sData.learning.velo.veloMeanMax);
sData.learning.velo.veloMeanMinPosBin = find(sData.behavior.binning.meanVeloBinned==sData.learning.velo.veloMeanMin);
if sData.learning.velo.veloMeanMinPosBin <= ceil(RewardZoneSize / sData.behavior.meta.binSize)
    sData.learning.velo.IsMeanMinInRewardZone = 'yes';
else
    sData.learning.velo.IsMeanMinInRewardZone = 'no';
end

sData.learning.velo.veloMax = max(sData.behavior.binning.veloBinned(1:nTrials,1:nBins),[],2);
sData.learning.velo.veloMin = min(sData.behavior.binning.veloBinned(1:nTrials,1:nBins),[],2);
sData.learning.velo.veloMaxPosBin = NaN(nTrials,1);
sData.learning.velo.veloMinPosBin = NaN(nTrials,1);
for i = 1:1:nTrials
    sData.learning.velo.veloMaxPosBin(i) = min(find(sData.behavior.binning.veloBinned(i,:) == sData.learning.velo.veloMax(i)));
    sData.learning.velo.veloMinPosBin(i) = min(find(sData.behavior.binning.veloBinned(i,:) == sData.learning.velo.veloMin(i)));
end
TrialsMinInRewardZone = find(sData.learning.velo.veloMinPosBin <= ceil(RewardZoneSize / sData.behavior.meta.binSize));
sData.learning.velo.MinInRewardZone = zeros(nTrials,1);
sData.learning.velo.MinInRewardZone(TrialsMinInRewardZone) = 1;
sData.learning.velo.MinInRewardZonePercent = length(TrialsMinInRewardZone)/nTrials*100;

% First and Second half od the recording VELO
sData.learning.velo.FirstSecondHalf = struct;
sData.learning.velo.FirstSecondHalf.MinBinPosFirst = mean(sData.learning.velo.veloMinPosBin(1:LapsAnalysed));
sData.learning.velo.FirstSecondHalf.MinBinPosLast = mean(sData.learning.velo.veloMinPosBin(nTrials-LapsAnalysed:nTrials));
sData.learning.velo.FirstSecondHalf.MinInRewardZonePercentFirst = sum(sData.learning.velo.MinInRewardZone(1:LapsAnalysed))/LapsAnalysed;
sData.learning.velo.FirstSecondHalf.MinInRewardZonePercentLast = sum(sData.learning.velo.MinInRewardZone(nTrials-LapsAnalysed:nTrials))/LapsAnalysed;
sData.learning.velo.FirstSecondHalf.Development_MinInRZPercentIncr = (sData.learning.velo.FirstSecondHalf.MinInRewardZonePercentLast - sData.learning.velo.FirstSecondHalf.MinInRewardZonePercentFirst)*100;
sData.learning.velo.FirstSecondHalf.Development_BinShiftFirstLast = sData.learning.velo.FirstSecondHalf.MinBinPosFirst - sData.learning.velo.FirstSecondHalf.MinBinPosLast;

% Save file to same path where LV files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% for copying data to excel
copy = NaN(1,9);
copy(1,1) = sData.learning.lick.LickPerformance_LateMiddle;
copy(1,2) = sData.learning.lick.LickPerformance_LateMiddleNorm;
copy(1,3) = sData.learning.lick.LickPerformance_LateEarly;
copy(1,4) = sData.learning.lick.LickPerformance_LateEarlyNorm;
copy(1,5) = sData.learning.lick.FirstSecondHalf.Development;
copy(1,7) = sData.learning.velo.MinInRewardZonePercent;
copy(1,8) = sData.learning.velo.FirstSecondHalf.Development_BinShiftFirstLast;
copy(1,9) = sData.learning.velo.FirstSecondHalf.Development_MinInRZPercentIncr;
copy(1,10) = sData.learning.velo.veloMeanMinPosBin;
sData.sessionInfo.sessionID

end