
AligmentCorr = struct;

nBin = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
roiStart = 1;
roiEnd = nROIs;
MotCorrShifts = sData.imdata.binned.MeanMotCorrBinned;
MeanRoiAct = NaN(nROIs,nBin);

% Mean ROI activity calc
for roi = roiStart:1:roiEnd
    MeanRoiAct(roi,1:nBin)= nanmean(sData.imdata.binned.RoidFF_SR{roi},1);
end

% normalization of MotCorr shifts and MeanRoiAct, and subtract MotCorrShifht from MeanRoiAct
MaxMotCorr = max(MotCorrShifts);
MinMotCorr = min(MotCorrShifts);
MotCorrShiftsNorm = (MotCorrShifts - MinMotCorr) ./  (MaxMotCorr - MinMotCorr); 

% normalization of MeanRoiAct
MaxMeanRoiAct = max(MeanRoiAct,[],2);
MinMeanRoiAct = min(MeanRoiAct,[],2);
MeanRoiActNorm =  (MeanRoiAct - MinMeanRoiAct) ./  (MaxMeanRoiAct - MinMeanRoiAct); 

MeanRoiActNormSubtr = NaN(nROIs,nBin);
NegModRois = [1 2 3 9 15 17];  %subtract MotCorrShiftsNorm
PosModRois = [4 5 6 7 8 11 12 14]; % add MotCorrShiftsNorm
% calc negModRois
for i=1:1:length(NegModRois)
    MeanRoiActNormSubtr(NegModRois(i),:) = MeanRoiActNorm(NegModRois(i),:) - MotCorrShiftsNorm(1,:);
end

% calc posModRois
for i=1:1:length(PosModRois)
    MeanRoiActNormSubtr(PosModRois(i),:) = MeanRoiActNorm(PosModRois(i),:) + MotCorrShiftsNorm(1,:);
end




% motion corretion vector error data was shifted with 2 bind backward to align with imaging data error
MeanRoiActNormSubtr2 = NaN(nROIs,nBin);
NegModRois = [1 2 3 9 15 17];  %subtract MotCorrShiftsNorm
PosModRois = [4 5 6 7 8 11 12 14]; % add MotCorrShiftsNorm
% calc negModRois
for i=1:1:length(NegModRois)
    MeanRoiActNormSubtr2(NegModRois(i),:) = MeanRoiActNorm(NegModRois(i),:) - MotCorrShiftNorm2binSh(1,:);
end

% calc posModRois
for i=1:1:length(PosModRois)
    MeanRoiActNormSubtr2(PosModRois(i),:) = MeanRoiActNorm(PosModRois(i),:) + MotCorrShiftNorm2binSh(1,:);
end
