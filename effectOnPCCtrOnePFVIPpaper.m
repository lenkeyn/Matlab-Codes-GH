function sData = effectOnPCCtrOnePFVIPpaper(sData)

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
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff.PlaceCellsWithOnePFList; % used those ROIs which are categorized as place cells in opto-off trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
    
% generate datastruct
sData.effectPCinCtrOnePFFinal = struct;
sData.effectPCinCtrOnePFFinal.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline)
sData.effectPCinCtrOnePFFinal.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCinCtrOnePFFinal.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCinCtrOnePFFinal.PFPeakPosBinOffOnAfter = NaN(nROIs,3);
sData.effectPCinCtrOnePFFinal.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3);
sData.effectPCinCtrOnePFFinal.meanActOutPF = NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCinCtrOnePFFinal.peakActInPFMeanOutRatio = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF

% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (in some ROIs there are more than one place fields), search for peak in position tuning curve
    
    %%% Opto-off
    PlaceFieldStartBinOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,1:20);
    PlaceFieldStartBinOff = PlaceFieldStartBinOffPre(~isnan(PlaceFieldStartBinOffPre));
    PlaceFieldLengthOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOff = PlaceFieldLengthOffPre(~isnan(PlaceFieldLengthOffPre));
    if sum(~isnan(PlaceFieldStartBinOffPre))== 0
       continue 
    end
    Peak = NaN(length(PlaceFieldLengthOff),1);
    for m = 1:1:length(PlaceFieldLengthOff) % if more then on place fields occurs
        Peak(m) = max(PosTuningOff(i,PlaceFieldStartBinOff(m):PlaceFieldStartBinOff(m)+PlaceFieldLengthOff(m)-1));
    end
    PeakOff = max(Peak);
    PeakStartPosOff = PlaceFieldStartBinOff(find(PeakOff));
    PeakLengthOff = PlaceFieldLengthOff(find(PeakOff));
    PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakOff));
    if PeakPosOff < DiscardBinsAtBegininng % if main peak is in the first bins which is to be discarded, discard this ROI from analysis
        continue
    end
    %Calculating effect size
    PosTuningForBaseline = PosTuningOff(i,1:nBins);
    PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
    PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field
    meanActOutPFOff = nanmean(PosTuningForBaseline);
    peakActInPFMeanOutRatioOff = PeakOff/meanActOutPFOff;
    
    %%% Opto-on
    PlaceFieldStartBinOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOn = PlaceFieldStartBinOnPre(~isnan(PlaceFieldStartBinOnPre));
    PlaceFieldLengthOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOn = PlaceFieldLengthOnPre(~isnan(PlaceFieldLengthOnPre));
    if sum(~isnan(PlaceFieldStartBinOnPre))== 0 % if there is no peak, take value from that position where peak was in Opto-off
        PeakOn = PosTuningOn(i,PeakPosOff);
        PeakStartPosOn = NaN;
        PeakLengthOn = NaN;
        PeakPosOn = NaN;
        %Calculating effect size
        PosTuningForBaseline = PosTuningOn(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field
        meanActOutPFOn = nanmean(PosTuningForBaseline);
        peakActInPFMeanOutRatioOn = PeakOn/meanActOutPFOn;
    else    
        Peak = NaN(length(PlaceFieldStartBinOn),1);
        for m = 1:1:length(PlaceFieldStartBinOn) % if more then on place fields occurs
            Peak(m) = max(PosTuningOn(i,PlaceFieldStartBinOn(m):PlaceFieldStartBinOn(m)+PlaceFieldLengthOn(m)-1));
        end
        PeakOn = max(Peak);
        PeakStartPosOn = PlaceFieldStartBinOn(find(PeakOn));
        PeakLengthOn = PlaceFieldStartBinOn(find(PeakOn));
        PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakOn));
        %Calculating effect size
        PosTuningForBaseline = PosTuningOn(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosOn:PeakStartPosOn+PeakLengthOn) = NaN; % discard data inside place field
        meanActOutPFOn = nanmean(PosTuningForBaseline);
        peakActInPFMeanOutRatioOn = PeakOn/meanActOutPFOn;
    end
    PeakCtrBinOn = PosTuningOn(i,PeakPosOff);
    
    %%% Opto-after
    PlaceFieldStartBinAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinAfter = PlaceFieldStartBinAfterPre(~isnan(PlaceFieldStartBinAfterPre));
    PlaceFieldLengthAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,:);
    PlaceFieldLengthAfter = PlaceFieldLengthAfterPre(~isnan(PlaceFieldLengthAfterPre));
    if sum(~isnan(PlaceFieldStartBinAfterPre))== 0 % if there is no peak, take value from that position where peak was in Opto-off
        PeakAfter = PosTuningAfter(i,PeakPosOff);
        PeakStartPosAfter = NaN;
        PeakLengthAfter = NaN;
        PeakPosAfter = NaN;
        %Calculating effect size
        PosTuningForBaseline = PosTuningAfter(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosOff:PeakStartPosOff+PeakLengthOff) = NaN; % discard data inside place field
        meanActOutPFAfter = nanmean(PosTuningForBaseline);
        peakActInPFMeanOutRatioAfter = PeakAfter/meanActOutPFAfter;
    else  
        Peak = NaN(length(PlaceFieldStartBinAfter),1);
        for m = 1:1:length(PlaceFieldStartBinAfter) % if more then on place fields occurs
            Peak(m) = max(PosTuningAfter(i,PlaceFieldStartBinAfter(m):PlaceFieldStartBinAfter(m)+PlaceFieldLengthAfter(m)-1));
        end
        PeakAfter = max(Peak);
        PeakStartPosAfter = PlaceFieldStartBinAfter(find(PeakAfter));
        PeakLengthAfter = PlaceFieldLengthAfter(find(max(Peak)));
        PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAfter));     
        %Calculating effect size
        PosTuningForBaseline = PosTuningAfter(i,1:nBins);
        PosTuningForBaseline(1:DiscardBinsAtBegininng) = NaN; % discard data in the beginning of trial
        PosTuningForBaseline(PeakStartPosAfter:PeakStartPosAfter+PeakLengthAfter) = NaN; % discard data inside place field
        meanActOutPFAfter = nanmean(PosTuningForBaseline);
        peakActInPFMeanOutRatioAfter = PeakAfter/meanActOutPFAfter;
    end
    PeakCtrBinAfter = PosTuningAfter(i,PeakPosOff);

    %%% save data for all protocol
    sData.effectPCinCtrOnePFFinal.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
    sData.effectPCinCtrOnePFFinal.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
    sData.effectPCinCtrOnePFFinal.PFPeakAmplOffOnAfter(i,1:3) = [PeakOff PeakOn PeakAfter];
    sData.effectPCinCtrOnePFFinal.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
    sData.effectPCinCtrOnePFFinal.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakOff PeakCtrBinOn PeakCtrBinAfter];
    sData.effectPCinCtrOnePFFinal.meanActOutPF(i,1:3) = [meanActOutPFOff meanActOutPFOn meanActOutPFAfter]; 
    sData.effectPCinCtrOnePFFinal.peakActInPFMeanOutRatio(i,1:3) = [peakActInPFMeanOutRatioOff peakActInPFMeanOutRatioOn peakActInPFMeanOutRatioAfter]; 
end

sData.effectPCinCtrOnePFFinal.note = 'Place cells detected in opto-off protocol';
sData.effectPCinCtrOnePFFinal.onlyValues = struct;
sData.effectPCinCtrOnePFFinal.onlyValues.PFStartPosBinOffOnAfter = sData.effectPCinCtrOnePFFinal.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
sData.effectPCinCtrOnePFFinal.onlyValues.PFLengthBinOffOnAfter = sData.effectPCinCtrOnePFFinal.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal.onlyValues.PFPeakAmplOffOnAfter = sData.effectPCinCtrOnePFFinal.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal.onlyValues.PFPeakPosBinOffOnAfter = sData.effectPCinCtrOnePFFinal.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal.PFPeakPosBinOffOnAfter(:,2))),1:3);
sData.effectPCinCtrOnePFFinal.onlyValues.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCinCtrOnePFFinal.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCinCtrOnePFFinal.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCinCtrOnePFFinal.onlyValues.meanActOutPF = sData.effectPCinCtrOnePFFinal.meanActOutPF(find(~isnan(sData.effectPCinCtrOnePFFinal.meanActOutPF(:,1))),1:3);
sData.effectPCinCtrOnePFFinal.onlyValues.peakActInPFMeanOutRatio = sData.effectPCinCtrOnePFFinal.peakActInPFMeanOutRatio(find(~isnan(sData.effectPCinCtrOnePFFinal.peakActInPFMeanOutRatio(:,1))),1:3);

% Save file to same path where other files can be found 
save(savePath,'sData');

end

