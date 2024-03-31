function sData = effectOnPCCtrVIPpaper2(sData)

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
DiscardBinsAtBegininng = ceil(sData.behavior.opto.optoStimStart) + 5; % there is no optical stimulation in the beginning of the trials, plus there is a delay (5 bins approx) even if laser is on , discard first bins data
PFPosJitter = 5; % bin, allowed jitter to consider PFs as same in different protocols
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells; % used those ROIs which are categorized as place cells in opto-off trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
    
% generate datastruct
sData.effectPCinCtrFinal2 = struct;
sData.effectPCinCtrFinal2.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3); % peak amplitude at that position which is the peak in Opto-off, calculated from pos tuning curve
sData.effectPCinCtrFinal2.meanActOutPFCtrPos = NaN(nROIs,3); % mean Ca activity outside of the PF, PF is determined in OPto-off
sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutCtrPos = NaN(nROIs,3); 
sData.effectPCinCtrFinal2.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline), can be different in different protocols
sData.effectPCinCtrFinal2.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCinCtrFinal2.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCinCtrFinal2.PFPeakPosBinOffOnAfter = NaN(nROIs,3); 
sData.effectPCinCtrFinal2.meanActOutPFNewPos= NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutNewPos = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF

% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0
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
    if sum((PeakStartPosOff - PFPosJitter < PotPlaceFieldStartBinAfter)&(PotPlaceFieldStartBinAfter < PeakStartPosOff + PFPosJitter))>0 % if there is a potential peak (PF starting) in opto on around the same location (+- 5 bin, PFPosJitter) as in opto off, use that peak
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
    % PF is determined in OPto-off
    sData.effectPCinCtrFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
    sData.effectPCinCtrFinal2.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
    sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
    % PF is allowed to shift +-5 bins in OPto-off and after protocols,
    % check values in these new PFs
    sData.effectPCinCtrFinal2.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
    sData.effectPCinCtrFinal2.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
    sData.effectPCinCtrFinal2.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
    sData.effectPCinCtrFinal2.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
    sData.effectPCinCtrFinal2.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
    sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
end

sData.effectPCinCtrFinal2.note = 'Place cells detected in opto-off protocol';
sData.effectPCinCtrFinal2.onlyValues = struct;
sData.effectPCinCtrFinal2.onlyValues.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCinCtrFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCinCtrFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCinCtrFinal2.onlyValues.meanActOutPFCtrPos = sData.effectPCinCtrFinal2.meanActOutPFCtrPos(find(~isnan(sData.effectPCinCtrFinal2.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCinCtrFinal2.onlyValues.RatioPeakActInPFMeanOutCtrPos = sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);

sData.effectPCinCtrFinal2.onlyValues.PFStartPosBinOffOnAfter = sData.effectPCinCtrFinal2.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrFinal2.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
sData.effectPCinCtrFinal2.onlyValues.PFLengthBinOffOnAfter = sData.effectPCinCtrFinal2.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCinCtrFinal2.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrFinal2.onlyValues.PFPeakAmplOffOnAfter = sData.effectPCinCtrFinal2.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCinCtrFinal2.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCinCtrFinal2.onlyValues.PFPeakPosBinOffOnAfter = sData.effectPCinCtrFinal2.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrFinal2.PFPeakPosBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrFinal2.onlyValues.meanActOutPFNewPos = sData.effectPCinCtrFinal2.meanActOutPFNewPos(find(~isnan(sData.effectPCinCtrFinal2.meanActOutPFNewPos(:,1))),1:3);
sData.effectPCinCtrFinal2.onlyValues.RatioPeakActInPFMeanOutNewPos = sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCinCtrFinal2.RatioPeakActInPFMeanOutNewPos(:,1))),1:3);

% Save file to same path where other files can be found 
save(savePath,'sData');

end


