function sData = effectOnPCLandmarkCellsVIPpaper3(sData)

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
%DiscardBinsAtBegininng = ceil(sData.behavior.opto.optoStimStart) + 5; % there is no optical stimulation in the beginning of the trials, plus there is a delay (5 bins approx) even if laser is on , discard first bins data
%PFPosJitter = 5; % bin, allowed jitter to consider PFs as same in different protocols

% position on Velcro oand HotGlue Spike landmarks:
LandmarkArea = floor(18/sData.behavior.meta.binSize); % 20 cm (changed from 4 cm), search for landmark cells where the landmark was presented and a little addition
VelcroStart1 = ceil(26/sData.behavior.meta.binSize); % 
VelcroEnd1 = VelcroStart1 + LandmarkArea;
VelcroStart2 = floor(86/sData.behavior.meta.binSize);
VelcroEnd2 = VelcroStart2 + LandmarkArea;
HGSStart1 = ceil(46/sData.behavior.meta.binSize);
HGSEnd1 = HGSStart1 + LandmarkArea;
HGSStart2 = ceil(66/sData.behavior.meta.binSize);
HGSEnd2 = HGSStart2 + LandmarkArea; % delay in Ca transient due to previous Ca transient evoked by first HGS

nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
LandmarkCell = sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.LandmarkCellinOptoOfforOptoOn; % used those ROIs which are categorized as place cells in opto-off or opto-on trials based on Mao criteria
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
    
% generate datastruct
sData.effectPCLandmarkCellFinal3 = struct;
sData.effectPCLandmarkCellFinal3.PFPeakAmplatCtrPeakPosOffOnAfter = NaN(nROIs,3); % peak amplitude at that position which is the peak in Opto-off, calculated from pos tuning curve
sData.effectPCLandmarkCellFinal3.meanActOutPFCtrPos = NaN(nROIs,3); % mean Ca activity outside of the PF, PF is determined in OPto-off
sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutCtrPos = NaN(nROIs,3); 
sData.effectPCLandmarkCellFinal3.PFStartPosBinOffOnAfter = NaN(nROIs,3);  % Starting bin of the PF (where pos tuning curve first goes above 1/3 of the peak compared to baseline), can be different in different protocols
sData.effectPCLandmarkCellFinal3.PFLengthBinOffOnAfter = NaN(nROIs,3); % length of PF (number of bins where pos tuning curve are above 1/3 of the peak compared to baseline)
sData.effectPCLandmarkCellFinal3.PFPeakAmplOffOnAfter = NaN(nROIs,3); % peak amplitude of PF , calculated from pos tuning curve
sData.effectPCLandmarkCellFinal3.PFPeakPosBinOffOnAfter = NaN(nROIs,3); 
sData.effectPCLandmarkCellFinal3.meanActOutPFNewPos= NaN(nROIs,3); % mean Ca activity outside of the PF (discard beginning of the track data, where no stimulation in opto-on trials)
sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutNewPos = NaN(nROIs,3); % ratio of peak activity inside/ mean activity outside the PF

% generate a common PlaceFieldStartBin, Length, peak amplitude, peaak position array for the three protocols (some cells are place cells in one, but not in another protocol.)
for i = 1:1:nROIs
    if sum(ismember(LandmarkCell,i))== 0
        continue
    end
    m = 0;
    % Which place field is the largest (in some ROIs there are more than one place fields), search for peak in position tuning curve
    if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.LandmarkCellsList,i))>0 % if it si a Landmark cell in opto-off
        m = 1;
        %%% Opto-off
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
            % Velcro cells must start to peak up between bin 14-16 and 44-46,
            % search the largest peak
            PeakAmplOff = max([PosTuningOff(i,VelcroStart1:VelcroEnd1) PosTuningOff(i,VelcroStart2:VelcroEnd2)]);
            PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff)); 
            PosTuningForBaseline = PosTuningOff(i,1:nBins);
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; % discard data in the trial until second velcro peak is done, in about 15 bins
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.HGSCellsList,i))>0
            % HGS cells must start to peak up  between bin 24-26 and 34-37
            PeakAmplOff = max([PosTuningOff(i,HGSStart1:HGSEnd1) PosTuningOff(i,HGSStart2:HGSEnd2)]);
            PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff));
            PosTuningForBaseline = PosTuningOff(i,1:nBins);
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; % discard data in the trial until second velcro peak is done
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
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.VelcroCellsList,i))>0 
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; 
        else
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; 
        end
        meanActOutPFCtrPosOn = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFCtrPosMeanOutOn = PeakAmplOnAtCtrPeakPos/meanActOutPFCtrPosOn;
        % check if the peak position has changed or not
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
            PeakAmplOn = max([PosTuningOn(i,VelcroStart1:VelcroEnd1) PosTuningOn(i,VelcroStart2:VelcroEnd2)]);
            PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn)); 
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.HGSCellsList,i))>0 % HGS cells 
            PeakAmplOn = max([PosTuningOn(i,HGSStart1:HGSEnd1) PosTuningOn(i,HGSStart2:HGSEnd2)]);
            PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn));
        end

        %%% Opto-after
        % take peak value from that position where peak was in Opto-off
        PeakAmplAfterAtCtrPeakPos = PosTuningAfter(i,PeakPosOff);
        %Calculating effect size
        PosTuningForBaseline = PosTuningAfter(i,1:nBins);
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.VelcroCellsList,i))>0 
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; 
        else
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; 
        end
        meanActOutPFCtrPosAfter = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFCtrPosMeanOutAfter = PeakAmplAfterAtCtrPeakPos/meanActOutPFCtrPosAfter;
        % check if the peak position has changed or not
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.VelcroCellsList,i))>0 % if it is a Velcro Cell
            PeakAmplAfter = max([PosTuningAfter(i,VelcroStart1:VelcroEnd1) PosTuningAfter(i,VelcroStart2:VelcroEnd2)]);
            PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter)); 
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOff.HGSCellsList,i))>0 % HGS cells 
            PeakAmplAfter = max([PosTuningAfter(i,HGSStart1:HGSEnd1) PosTuningAfter(i,HGSStart2:HGSEnd2)]);
            PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter));
        end
    %%% if it si a Landmark cell in opto-off    
    elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.LandmarkCellsList,i))>0 % if it si a Landmark cell in opto-off
        m = 2;
        %%% Opto-on
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.VelcroCellsList,i))>0 % if it is a Velcro Cell
            % Velcro cells must start to peak up between bin 14-16 and 44-46,
            % search the largest peak
            PeakAmplOn = max([PosTuningOn(i,VelcroStart1:VelcroEnd1) PosTuningOn(i,VelcroStart2:VelcroEnd2)]);
            PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn)); 
            PosTuningForBaseline = PosTuningOn(i,1:nBins);
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; % discard data in the trial until second velcro peak is done, in about 15 bins
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.HGSCellsList,i))>0
            % HGS cells must start to peak up  between bin 24-26 and 34-37
            PeakAmplOn = max([PosTuningOn(i,HGSStart1:HGSEnd1) PosTuningOn(i,HGSStart2:HGSEnd2)]);
            PeakPosOn = min(find(PosTuningOn(i,1:nBins)==PeakAmplOn));
            PosTuningForBaseline = PosTuningOn(i,1:nBins);
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; % discard data in the trial until second velcro peak is done
        end
        %Calculating effect size
        meanActOutPFOn = nanmean(PosTuningForBaseline);
        peakActInPFMeanOutRatioOn = PeakAmplOn/meanActOutPFOn;
    
        %%% Opto-off
        % take peak value from that position where peak was in Opto-off
        PeakAmplOffAtCtrPeakPos = PosTuningOff(i,PeakPosOn);
        %Calculating effect size
        PosTuningForBaseline = PosTuningOff(i,1:nBins);
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.VelcroCellsList,i))>0 
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; 
        else
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; 
        end
        meanActOutPFCtrPosOff = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFCtrPosMeanOutOff = PeakAmplOffAtCtrPeakPos/meanActOutPFCtrPosOff;
        % check if the peak position has changed or not
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.VelcroCellsList,i))>0 % if it is a Velcro Cell
            PeakAmplOff = max([PosTuningOff(i,VelcroStart1:VelcroEnd1) PosTuningOff(i,VelcroStart2:VelcroEnd2)]);
            PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff)); 
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.HGSCellsList,i))>0 % HGS cells 
            PeakAmplOff = max([PosTuningOff(i,HGSStart1:HGSEnd1) PosTuningOff(i,HGSStart2:HGSEnd2)]);
            PeakPosOff = min(find(PosTuningOff(i,1:nBins)==PeakAmplOff));
        end

        %%% Opto-after
        % take peak value from that position where peak was in Opto-on
        PeakAmplAfterAtCtrPeakPos = PosTuningAfter(i,PeakPosOn);
        %Calculating effect size
        PosTuningForBaseline = PosTuningAfter(i,1:nBins);
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.VelcroCellsList,i))>0 
            PosTuningForBaseline(1:VelcroEnd2 + 15) = NaN; 
        else
            PosTuningForBaseline(1:HGSEnd2 + 15) = NaN; 
        end
        meanActOutPFCtrPosAfter = nanmean(PosTuningForBaseline);
        RATIOpeakActInPFCtrPosMeanOutAfter = PeakAmplAfterAtCtrPeakPos/meanActOutPFCtrPosAfter;
        % check if the peak position has changed or not
        if sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.VelcroCellsList,i))>0 % if it is a Velcro Cell
            PeakAmplAfter = max([PosTuningAfter(i,VelcroStart1:VelcroEnd1) PosTuningAfter(i,VelcroStart2:VelcroEnd2)]);
            PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter)); 
        elseif  sum(ismember(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.OptoOn.HGSCellsList,i))>0 % HGS cells 
            PeakAmplAfter = max([PosTuningAfter(i,HGSStart1:HGSEnd1) PosTuningAfter(i,HGSStart2:HGSEnd2)]);
            PeakPosAfter = min(find(PosTuningAfter(i,1:nBins)==PeakAmplAfter));
        end
    end
    
 %%% save data for all protocol
    % PF is determined in OPto-off
    if m == 1
        sData.effectPCLandmarkCellFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOnAtCtrPeakPos PeakAmplAfterAtCtrPeakPos];
        sData.effectPCLandmarkCellFinal3.meanActOutPFCtrPos(i,1:3) = [meanActOutPFOff meanActOutPFCtrPosOn meanActOutPFCtrPosAfter];
        sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFCtrPosMeanOutOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
    elseif m == 2
        sData.effectPCLandmarkCellFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(i,1:3) = [PeakAmplOffAtCtrPeakPos PeakAmplOn PeakAmplAfterAtCtrPeakPos];
        sData.effectPCLandmarkCellFinal3.meanActOutPFCtrPos(i,1:3) = [meanActOutPFCtrPosOff meanActOutPFOn meanActOutPFCtrPosAfter];
        sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutCtrPos(i,1:3) = [RATIOpeakActInPFCtrPosMeanOutOff peakActInPFMeanOutRatioOn RATIOpeakActInPFCtrPosMeanOutAfter]; 
    end
    % PF is allowed to shift +-5 bins in OPto-off and after protocols,
    % check values in these new PFs
    %sData.effectPCLandmarkCellFinal3.PFStartPosBinOffOnAfter(i,1:3) = [PeakStartPosOff PeakStartPosOn PeakStartPosAfter]; 
    %sData.effectPCLandmarkCellFinal3.PFLengthBinOffOnAfter(i,1:3) = [PeakLengthOff PeakLengthOn PeakLengthAfter];
    sData.effectPCLandmarkCellFinal3.PFPeakAmplOffOnAfter(i,1:3) = [PeakAmplOff PeakAmplOn PeakAmplAfter];
    sData.effectPCLandmarkCellFinal3.PFPeakPosBinOffOnAfter(i,1:3) = [PeakPosOff PeakPosOn PeakPosAfter]; 
    %sData.effectPCLandmarkCellFinal3.meanActOutPFNewPos(i,1:3) = [meanActOutPFOff meanActOutPFNewPosOn meanActOutPFNewPosAfter]; 
    %sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutNewPos(i,1:3) = [peakActInPFMeanOutRatioOff RATIOpeakActInPFNewPosMeanOutOn RATIOpeakActInPFNewPosMeanOutAfter]; 
    
end

sData.effectPCLandmarkCellFinal3.note = 'Landmark cells were detected in opto-off protocol';
sData.effectPCLandmarkCellFinal3.onlyValues = struct;
sData.effectPCLandmarkCellFinal3.onlyValues.PFPeakAmplatCtrPeakPosOffOnAfter = sData.effectPCLandmarkCellFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal3.PFPeakAmplatCtrPeakPosOffOnAfter(:,1))),1:3);
sData.effectPCLandmarkCellFinal3.onlyValues.meanActOutPFCtrPos = sData.effectPCLandmarkCellFinal3.meanActOutPFCtrPos(find(~isnan(sData.effectPCLandmarkCellFinal3.meanActOutPFCtrPos(:,1))),1:3);
sData.effectPCLandmarkCellFinal3.onlyValues.RatioPeakActInPFMeanOutCtrPos = sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutCtrPos(find(~isnan(sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutCtrPos(:,1))),1:3);

%sData.effectPCLandmarkCellFinal3.onlyValues.PFStartPosBinOffOnAfter = sData.effectPCLandmarkCellFinal3.PFStartPosBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal3.PFStartPosBinOffOnAfter(:,2))),1:3); % value should not be NaN in Opto-On trials
%sData.effectPCLandmarkCellFinal3.onlyValues.PFLengthBinOffOnAfter = sData.effectPCLandmarkCellFinal3.PFLengthBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal3.PFLengthBinOffOnAfter(:,2))),1:3);
sData.effectPCLandmarkCellFinal3.onlyValues.PFPeakAmplOffOnAfter = sData.effectPCLandmarkCellFinal3.PFPeakAmplOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal3.PFPeakAmplOffOnAfter(:,2))),1:3);
sData.effectPCLandmarkCellFinal3.onlyValues.PFPeakPosBinOffOnAfter = sData.effectPCLandmarkCellFinal3.PFPeakPosBinOffOnAfter(find(~isnan(sData.effectPCLandmarkCellFinal3.PFPeakPosBinOffOnAfter(:,2))),1:3);
%sData.effectPCLandmarkCellFinal3.onlyValues.meanActOutPFNewPos = sData.effectPCLandmarkCellFinal3.meanActOutPFNewPos(find(~isnan(sData.effectPCLandmarkCellFinal3.meanActOutPFNewPos(:,2))),1:3);
%sData.effectPCLandmarkCellFinal3.onlyValues.RatioPeakActInPFMeanOutNewPos = sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutNewPos(find(~isnan(sData.effectPCLandmarkCellFinal3.RatioPeakActInPFMeanOutNewPos(:,2))),1:3);


% Save file to same path where other files can be found 
save(savePath,'sData');

end


