function sData = effectOnPCLandmarkCellsVIPpaper2(sData)

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
% position on Velcr oand HotGlue Spike landmarks:
Velcro1Start = floor(28/sData.behavior.meta.binSize); % 28 cm
Velcro1End = floor(44/sData.behavior.meta.binSize); % 44 cm, peak should be reached until then
Velcro2Start = floor(88/sData.behavior.meta.binSize); % 88 cm
Velcro2End = floor(104/sData.behavior.meta.binSize); % 104 cm, peak should be reached until then
HGS1Start = floor(48/sData.behavior.meta.binSize); % 48 cm
HGS1End = floor(64/sData.behavior.meta.binSize); % 64 cm, peak should be reached until then
HGS2Start = floor(68/sData.behavior.meta.binSize); % 68 cm
HGS2End = floor(84/sData.behavior.meta.binSize); % 84 cm, peak should be reached until then
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.LandmarkCellsList; % used those ROIs which are categorized as place cells in opto-off trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
    
% generate datastruct
sData.effectPCLandmarkCellFinal2 = struct;
sData.effectPCLandmarkCellFinal2.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3); % peak amplitude at that position which is the peak in Opto-off, calculated from pos tuning curve
sData.effectPCLandmarkCellFinal2.meanActOutPFCtrPos = NaN(nROIs,3); % mean Ca activity outside of the PF, PF is determined in OPto-off
sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutCtrPos = NaN(nROIs,3); 
sData.effectPCLandmarkCellFinal2.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline), can be different in different protocols
sData.effectPCLandmarkCellFinal2.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCLandmarkCellFinal2.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCLandmarkCellFinal2.PFPeakPosBinOffOnAfter = NaN(nROIs,3); 
sData.effectPCLandmarkCellFinal2.meanActOutPFNewPos= NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutNewPos = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF

% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (in some ROIs there are more than one place fields), search for peak in position tuning curve
    
    %%% Opto-off
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
        % Velcro cells must start to peak up between bin 14-16 and 44-46,
        % search the largest peak
        PeakAmplOff = max([PosTuningOff(i,Velcro1Start:Velcro1End) PosTuningOff(i,Velcro2Start:Velcro2End)]);
        PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff)); 
        PosTuningForBaseline = PosTuningOff(i,1:nBins);
        PosTuningForBaseline(1:Velcro2End + 15) = NaN; % discard data in the trial until second velcro peak is done, in about 15 bins
    elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.HGSCellsList,i))>0
        % HGS cells must start to peak up  between bin 24-26 and 34-37
        PeakAmplOff = max([PosTuningOff(i,HGS1Start:HGS1End) PosTuningOff(i,HGS2Start:HGS2End)]);
        PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff));
        PosTuningForBaseline = PosTuningOff(i,1:nBins);
        PosTuningForBaseline(1:HGS2End + 15) = NaN; % discard data in the trial until second velcro peak is done
    end
    %Calculating effect size
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
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.VelcroCellsList,i))>0 
        PosTuningForBaseline(1:Velcro2End + 15) = NaN; 
    else
        PosTuningForBaseline(1:HGS2End + 15) = NaN; 
    end
    meanActOutPFCtrPosOn = nanmean(PosTuningForBaseline);
    RATIOpeakActInPFCtrPosMeanOutOn = PeakAmplOnAtCtrPeakPos/meanActOutPFCtrPosOn;
    % check if the peak position has changed or not
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
        PeakAmplOn = max([PosTuningOn(i,Velcro1Start:Velcro1End) PosTuningOn(i,Velcro2Start:Velcro2End)]);
        PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn)); 
    elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.HGSCellsList,i))>0 % HGS cells 
        PeakAmplOn = max([PosTuningOn(i,HGS1Start:HGS1End) PosTuningOn(i,HGS2Start:HGS2End)]);
        PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn));
    end
            
    %%% Opto-after
    % take peak value from that position where peak was in Opto-off
    PeakAmplAfterAtCtrPeakPos = PosTuningAfter(i,PeakPosOff);
    %Calculating effect size
    PosTuningForBaseline = PosTuningAfter(i,1:nBins);
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.VelcroCellsList,i))>0 
        PosTuningForBaseline(1:Velcro2End + 15) = NaN; 
    else
        PosTuningForBaseline(1:HGS2End + 15) = NaN; 
    end
    meanActOutPFCtrPosAfter = nanmean(PosTuningForBaseline);
    RATIOpeakActInPFCtrPosMeanOutAfter = PeakAmplAfterAtCtrPeakPos/meanActOutPFCtrPosAfter;
    % check if the peak position has changed or not
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
        PeakAmplAfter = max([PosTuningAfter(i,Velcro1Start:Velcro1End) PosTuningAfter(i,Velcro2Start:Velcro2End)]);
        PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter)); 
    elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.HGSCellsList,i))>0 % HGS cells 
        PeakAmplAfter = max([PosTuningAfter(i,HGS1Start:HGS1End) PosTuningAfter(i,HGS2Start:HGS2End)]);
        PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter));
    end
     
    %%% save data for all protocol
    % PF is determined in OPto-off
    sData.effectPCLandmarkCellFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
    sData.effectPCLandmarkCellFinal2.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
    sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
    
    % PF is allowed to shift +-5 bins in OPto-off and after protocols,
    % check values in these new PFs
    %sData.effectPCLandmarkCellFinal2.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
    %sData.effectPCLandmarkCellFinal2.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
    sData.effectPCLandmarkCellFinal2.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
    sData.effectPCLandmarkCellFinal2.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
    %sData.effectPCLandmarkCellFinal2.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
    %sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
    
end

sData.effectPCLandmarkCellFinal2.note = 'Landmark cells were detected in opto-off protocol';
sData.effectPCLandmarkCellFinal2.onlyValues = struct;
sData.effectPCLandmarkCellFinal2.onlyValues.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCLandmarkCellFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal2.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCLandmarkCellFinal2.onlyValues.meanActOutPFCtrPos = sData.effectPCLandmarkCellFinal2.meanActOutPFCtrPos(find(~isnan(sData.effectPCLandmarkCellFinal2.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCLandmarkCellFinal2.onlyValues.RatioPeakActInPFMeanOutCtrPos = sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);

%sData.effectPCLandmarkCellFinal2.onlyValues.PFStartPosBinOffOnAfter = sData.effectPCLandmarkCellFinal2.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal2.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
%sData.effectPCLandmarkCellFinal2.onlyValues.PFLengthBinOffOnAfter = sData.effectPCLandmarkCellFinal2.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal2.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCLandmarkCellFinal2.onlyValues.PFPeakAmplOffOnAfter = sData.effectPCLandmarkCellFinal2.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal2.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCLandmarkCellFinal2.onlyValues.PFPeakPosBinOffOnAfter = sData.effectPCLandmarkCellFinal2.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal2.PFPeakPosBinOffOnAfter(:,2))),1:3);
%sData.effectPCLandmarkCellFinal2.onlyValues.meanActOutPFNewPos = sData.effectPCLandmarkCellFinal2.meanActOutPFNewPos(find(~isnan(sData.effectPCLandmarkCellFinal2.meanActOutPFNewPos(:,2))),1:3);
%sData.effectPCLandmarkCellFinal2.onlyValues.RatioPeakActInPFMeanOutNewPos = sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCLandmarkCellFinal2.RatioPeakActInPFMeanOutNewPos(:,2))),1:3);


% Save file to same path where other files can be found 
save(savePath,'sData');

end


