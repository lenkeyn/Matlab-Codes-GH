function sData = effectOnPCBeforeRewardVIPpaper3(sData)

% First you have to run effectOnPCCtrOnePFVIPpaper2

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects all place cells which are categorized as a place cell in control.
% compare peak amlitude and position of peak in control and opto trials
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols

%Set saving:
Folder = 'C:\MATLAB\SAVE\FinalEffect';
savePath = fullfile(Folder,sData.sessionInfo.fileID);
if ~isfolder(savePath) 
    mkdir(savePath)
end

%Set paramteres: 
DiscardBinsAtBegininng = round(30/sData.behavior.meta.binSize); % don't consider PCs which has peak in earlier location (speed and PM VIP activity goes up until 30 cm, and starts to decresase at 127 cm, 30 cm before reward)
RewardBinsStart = round(127/sData.behavior.meta.binSize);
PFPosJitter = 5; % bin, allowed jitter to consider PFs as same in different protocols
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff.PlaceCellsWithOnePFList; % used those ROIs which are categorized as place cells in opto-off trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
    
% generate datastruct
sData.effectPCvsRewardCellFinal2 = struct;
% PC which peak are between 30 cm and 127 cm:
sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3); % peak amplitude at that position which is the peak in Opto-off, calculated from pos tuning curve
sData.effectPCvsRewardCellFinal2.PC.meanActOutPFCtrPos = NaN(nROIs,3); % mean Ca activity outside of the PF, PF is determined in OPto-off
sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutCtrPos = NaN(nROIs,3); 
sData.effectPCvsRewardCellFinal2.PC.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline), can be different in different protocols
sData.effectPCvsRewardCellFinal2.PC.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCvsRewardCellFinal2.PC.PFPeakPosBinOffOnAfter = NaN(nROIs,3); 
sData.effectPCvsRewardCellFinal2.PC.meanActOutPFNewPos= NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutNewPos = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF
% PC which peak are between 127 cm and 157 cm:
sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3); % peak amplitude at that position which is the peak in Opto-off, calculated from pos tuning curve
sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFCtrPos = NaN(nROIs,3); % mean Ca activity outside of the PF, PF is determined in OPto-off
sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutCtrPos = NaN(nROIs,3); 
sData.effectPCvsRewardCellFinal2.RewardCell.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline), can be different in different protocols
sData.effectPCvsRewardCellFinal2.RewardCell.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakPosBinOffOnAfter = NaN(nROIs,3); 
sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFNewPos= NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutNewPos = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF

% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0  ||  sum(ismember(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (in some ROIs there are more than one place fields), search for peak in position tuning curve
    
    %%% Opto-off
    PlaceFieldStartBinOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOff = PlaceFieldStartBinOffPre(~isnan(PlaceFieldStartBinOffPre));
    PlaceFieldLengthOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOff = PlaceFieldLengthOffPre(~isnan(PlaceFieldLengthOffPre));
    Peak = NaN(length(PlaceFieldLengthOff),1);
    for m = 1:1:length(PlaceFieldLengthOff) % if more then on place fields occurs
        Peak(m) = max(PosTuningOff(i,PlaceFieldStartBinOff(m):PlaceFieldStartBinOff(m)+PlaceFieldLengthOff(m)-1));
    end
    PeakAmplOff = max(Peak);
    PeakStartPosOff = PlaceFieldStartBinOff(find(PeakAmplOff));
    PeakLengthOff = PlaceFieldLengthOff(find(PeakAmplOff));
    PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff));
    if PeakPosOff < DiscardBinsAtBegininng % if main peak is in the first bins which is to be discarded, discard this ROI from analysis
        continue
    end
    %Calculating effect size
    PosTuningForBaseline = PosTuningOff(i,1:nBins);
    PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
    PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field
    meanActOutPFOff = nanmean(PosTuningForBaseline);
    peakActInPFMeanOutRatioOff = PeakAmplOff/meanActOutPFOff;
    
    %%% Opto-on
    % previous code searched for the place fields in each protocol, but
    % later I just checked the change of the already detected pf in
    % control, because in some ROIs the red laser switch on caused that
    % large activation in the beginning that that was detected as the main place field.

    % take peak value from that position where peak was in Opto-off
    PeakAmplOnAtCtrPeakPos = PosTuningOn(i,PeakPosOff);
    %Calculating effect size
    PosTuningForBaseline = PosTuningOn(i,1:nBins);
    PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
    PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field based on opto-off
    meanActOutPFCtrPosOn = nanmean(PosTuningForBaseline);
    RATIOpeakActInPFCtrPosMeanOutOn = PeakAmplOnAtCtrPeakPos/meanActOutPFCtrPosOn;
    
    % search the peak value and its location in Opto-on within PF location in Opto-off (position can be little bit different as in Opto-off)  
    PotPlaceFieldStartBinOn = sData.imdata.MaoPC_Opto_dff.OptoOn.PotPlaceFieldStartBin(i,find(~isnan(sData.imdata.MaoPC_Opto_dff.OptoOn.PotPlaceFieldStartBin(i,:))));
    PotPlaceFieldLengthOn = sData.imdata.MaoPC_Opto_dff.OptoOn.PotPlaceFieldLength(i,find(~isnan(sData.imdata.MaoPC_Opto_dff.OptoOn.PotPlaceFieldLength(i,:))));
    if sum((PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinOn) & (PotPlaceFieldStartBinOn < PeakStartPosOff + PFPosJitter))>0 % if there is a potential peak (PF starting) in opto on around the same location (+- 5 bin, PFPosJitter) as in opto off, use that peak
        PlaceFieldStartBinOn = PotPlaceFieldStartBinOn(intersect(find(PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinOn),find(PeakStartPosOff + PFPosJitter > PotPlaceFieldStartBinOn)));    
        PlaceFieldLengthOn = PotPlaceFieldLengthOn(intersect(find(PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinOn),find(PeakStartPosOff + PFPosJitter > PotPlaceFieldStartBinOn)));    
        Peak = NaN(length(PlaceFieldStartBinOn),1);
        for m = 1:1:length(PlaceFieldStartBinOn) % if more then on place fields occurs
            Peak(m) = max(PosTuningOn(i,PlaceFieldStartBinOn(m):PlaceFieldStartBinOn(m)+PlaceFieldLengthOn(m)));
        end
        PeakAmplOn = max(Peak);
        PeakStartPosOn = PlaceFieldStartBinOn(find(PeakAmplOn));
        PeakLengthOn = PlaceFieldStartBinOn(find(PeakAmplOn));
        PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn));
        %Calculating effect size
        PosTuningForBaseline = PosTuningOn(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosOn:PeakStartPosOn+PeakLengthOn) = NaN; % discard data inside place field
        meanActOutPFNewPosOn = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFNewPosMeanOutOn = PeakAmplOn/meanActOutPFNewPosOn;
    else
        PeakAmplOn = NaN;
        PeakStartPosOn = NaN;
        PeakLengthOn = NaN;
        PeakPosOn = NaN;
        meanActOutPFNewPosOn = NaN;
        RATIOpeakActInPFNewPosMeanOutOn = NaN;
    end
        
    %%% Opto-after
    % take peak value from that position where peak was in Opto-off
    PeakAmplAfterAtCtrPeakPos = PosTuningAfter(i,PeakPosOff);
    %Calculating effect size
    PosTuningForBaseline = PosTuningAfter(i,1:nBins);
    PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
    PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field based on opto-off
    meanActOutPFCtrPosAfter = nanmean(PosTuningForBaseline);
    RATIOpeakActInPFCtrPosMeanOutAfter = PeakAmplAfterAtCtrPeakPos/meanActOutPFCtrPosAfter;
    
    % search the peak value and its location in Opto-after within PF location in Opto-off (position can be little bit different as in Opto-off)  
    PotPlaceFieldStartBinAfter = sData.imdata.MaoPC_Opto_dff.OptoAfter.PotPlaceFieldStartBin(i,find(~isnan(sData.imdata.MaoPC_Opto_dff.OptoAfter.PotPlaceFieldStartBin(i,:))));
    PotPlaceFieldLengthAfter = sData.imdata.MaoPC_Opto_dff.OptoAfter.PotPlaceFieldLength(i,find(~isnan(sData.imdata.MaoPC_Opto_dff.OptoAfter.PotPlaceFieldLength(i,:))));
    if sum((PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinAfter) & (PotPlaceFieldStartBinAfter < PeakStartPosOff + PFPosJitter))>0 % if there is a potential peak (PF starting) in opto on around the same location (+- 5 bin, PFPosJitter) as in opto off, use that peak
        PlaceFieldStartBinAfter = PotPlaceFieldStartBinAfter(intersect(find(PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinAfter),find(PeakStartPosOff + PFPosJitter > PotPlaceFieldStartBinAfter)));    
        PlaceFieldLengthAfter = PotPlaceFieldLengthAfter(intersect(find(PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinAfter),find(PeakStartPosOff + PFPosJitter > PotPlaceFieldStartBinAfter)));    
        Peak = NaN(length(PlaceFieldStartBinAfter),1);
        for m = 1:1:length(PlaceFieldStartBinAfter) % if more then on place fields occurs
            Peak(m) = max(PosTuningAfter(i,PlaceFieldStartBinAfter(m):PlaceFieldStartBinAfter(m)+PlaceFieldLengthAfter(m)));
        end
        PeakAmplAfter = max(Peak);
        PeakStartPosAfter = PlaceFieldStartBinAfter(find(PeakAmplAfter));
        PeakLengthAfter = PlaceFieldStartBinAfter(find(PeakAmplAfter));
        PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter));
        %Calculating effect size
        PosTuningForBaseline = PosTuningAfter(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosAfter:PeakStartPosAfter+PeakLengthAfter) = NaN; % discard data inside place field
        meanActOutPFNewPosAfter = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFNewPosMeanOutAfter = PeakAmplAfter/meanActOutPFNewPosAfter;
    else
        PeakAmplAfter = NaN;
        PeakStartPosAfter = NaN;
        PeakLengthAfter = NaN;
        PeakPosAfter = NaN;
        meanActOutPFNewPosAfter = NaN;
        RATIOpeakActInPFNewPosMeanOutAfter = NaN;
    end
    
    
    %%% save data for all protocol
    % if it is a normal PC:
    if PeakPosOff < RewardBinsStart
        % PF is determined in OPto-off
        sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
        sData.effectPCvsRewardCellFinal2.PC.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
        sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
        % PF is allowed to shift +-5 bins in OPto-off and after protocols,
        % check values in these new PFs
        sData.effectPCvsRewardCellFinal2.PC.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
        sData.effectPCvsRewardCellFinal2.PC.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
        sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
        sData.effectPCvsRewardCellFinal2.PC.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
        sData.effectPCvsRewardCellFinal2.PC.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
        sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
    else % if it is a reward cell
        % PF is determined in OPto-off
        sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
        sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
        sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
        % PF is allowed to shift +-5 bins in OPto-off and after protocols,
        % check values in these new PFs
        sData.effectPCvsRewardCellFinal2.RewardCell.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
        sData.effectPCvsRewardCellFinal2.RewardCell.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
        sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
        sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
        sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
        sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
    end
end

sData.effectPCvsRewardCellFinal2.note = 'Place cells detected in opto-off protocol';
sData.effectPCvsRewardCellFinal2.onlyValues = struct;

sData.effectPCvsRewardCellFinal2.onlyValues.PC.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.meanActOutPFCtrPos = sData.effectPCvsRewardCellFinal2.PC.meanActOutPFCtrPos(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.RatioPeakActInPFMeanOutCtrPos = sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.PFStartPosBinOffOnAfter = sData.effectPCvsRewardCellFinal2.PC.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
sData.effectPCvsRewardCellFinal2.onlyValues.PC.PFLengthBinOffOnAfter = sData.effectPCvsRewardCellFinal2.PC.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.PFPeakAmplOffOnAfter = sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.PFPeakPosBinOffOnAfter = sData.effectPCvsRewardCellFinal2.PC.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.PFPeakPosBinOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.meanActOutPFNewPos = sData.effectPCvsRewardCellFinal2.PC.meanActOutPFNewPos(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.meanActOutPFNewPos(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.PC.RatioPeakActInPFMeanOutNewPos = sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCvsRewardCellFinal2.PC.RatioPeakActInPFMeanOutNewPos(:,2))),1:3);

sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.meanActOutPFCtrPos = sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFCtrPos(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.RatioPeakActInPFMeanOutCtrPos = sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.PFStartPosBinOffOnAfter = sData.effectPCvsRewardCellFinal2.RewardCell.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.PFLengthBinOffOnAfter = sData.effectPCvsRewardCellFinal2.RewardCell.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.PFPeakAmplOffOnAfter = sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.PFPeakPosBinOffOnAfter = sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.PFPeakPosBinOffOnAfter(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.meanActOutPFNewPos = sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFNewPos(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.meanActOutPFNewPos(:,2))),1:3);
sData.effectPCvsRewardCellFinal2.onlyValues.RewardCell.RatioPeakActInPFMeanOutNewPos = sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCvsRewardCellFinal2.RewardCell.RatioPeakActInPFMeanOutNewPos(:,2))),1:3);

% Save file to same path where other files can be found 
save(savePath,'sData');

end


