function sData = effectOnPCCtrOnePFVIPpaper3(sData)

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
PlaceCells = sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.PlaceCellsWithOnePFinOptoOfforOn; % used those ROIs which are categorized as place cells in opto-off trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
   
% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0  %||  sum(ismember(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (in some ROIs there are more than one place fields), search for peak in position tuning curve
    
    %%% Opto-off
    PlaceFieldStartBinOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,:);
    if sum(~isnan(PlaceFieldStartBinOffPre))== 0
        continue
    else
        PlaceFieldStartBinOff = PlaceFieldStartBinOffPre(~isnan(PlaceFieldStartBinOffPre));
        PlaceFieldLengthOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,:);
        PlaceFieldLengthOff = PlaceFieldLengthOffPre(~isnan(PlaceFieldLengthOffPre));
        Peak = NaN(length(PlaceFieldLengthOff),1);
        for m = 1:1:length(PlaceFieldLengthOff) % if more then on place fields occurs
            Peak(m) = max(PosTuningOff(i,PlaceFieldStartBinOff(m):PlaceFieldStartBinOff(m)+PlaceFieldLengthOff(m)-1));
        end
        PeakAmplOff = max(Peak);
        PeakStartPosOff = PlaceFieldStartBinOff(find(PeakAmplOff==Peak));
        PeakLengthOff = PlaceFieldLengthOff(find(PeakAmplOff==Peak));
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
    end
    
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
        PeakStartPosOn = PlaceFieldStartBinOn(find(PeakAmplOn==Peak));
        PeakLengthOn = PlaceFieldStartBinOn(find(PeakAmplOn==Peak));
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
        PeakStartPosAfter = PlaceFieldStartBinAfter(find(PeakAmplAfter==Peak));
        PeakLengthAfter = PlaceFieldStartBinAfter(find(PeakAmplAfter==Peak));
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
    sData.effectPCinCtrOnePFFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
    sData.effectPCinCtrOnePFFinal3.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
    sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
    % PF is allowed to shift +-5 bins in OPto-off and after protocols,
    % check values in these new PFs
    sData.effectPCinCtrOnePFFinal3.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
    sData.effectPCinCtrOnePFFinal3.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
    sData.effectPCinCtrOnePFFinal3.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
    sData.effectPCinCtrOnePFFinal3.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
    sData.effectPCinCtrOnePFFinal3.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
    sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
end

sData.effectPCinCtrOnePFFinal3.note = 'Place cells detected in opto-off protocol';
sData.effectPCinCtrOnePFFinal3.onlyValues = struct;
sData.effectPCinCtrOnePFFinal3.onlyValues.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCinCtrOnePFFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.meanActOutPFCtrPos = sData.effectPCinCtrOnePFFinal3.meanActOutPFCtrPos(find(~isnan(sData.effectPCinCtrOnePFFinal3.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.RatioPeakActInPFMeanOutCtrPos = sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);

sData.effectPCinCtrOnePFFinal3.onlyValues.PFStartPosBinOffOnAfter = sData.effectPCinCtrOnePFFinal3.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal3.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
sData.effectPCinCtrOnePFFinal3.onlyValues.PFLengthBinOffOnAfter = sData.effectPCinCtrOnePFFinal3.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal3.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.PFPeakAmplOffOnAfter = sData.effectPCinCtrOnePFFinal3.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal3.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.PFPeakPosBinOffOnAfter = sData.effectPCinCtrOnePFFinal3.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal3.PFPeakPosBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.meanActOutPFNewPos = sData.effectPCinCtrOnePFFinal3.meanActOutPFNewPos(find(~isnan(sData.effectPCinCtrOnePFFinal3.meanActOutPFNewPos(:,2))),1:3);
sData.effectPCinCtrOnePFFinal3.onlyValues.RatioPeakActInPFMeanOutNewPos = sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCinCtrOnePFFinal3.RatioPeakActInPFMeanOutNewPos(:,2))),1:3);

% Save file to same path where other files can be found 
save(savePath,'sData');

end


