function sData = gainModulationInOutRatio(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.

% create four arrays of trials: opto-off-1, opto-off-2 (to calculate a control between two opto-off), opto-on, after-opto
%{
RandIndices = randsample(1:length(sData.behavior.opto.OptoOffTrialsIndices),floor(length(sData.behavior.opto.OptoOffTrialsIndices)/2));
OptoOffTrials1 = sort(sData.behavior.opto.OptoOffTrialsIndices(RandIndices));
OptoOffTrialsTemp =  sData.behavior.opto.OptoOffTrialsIndices;
OptoOffTrialsTemp(RandIndices) = NaN;
OptoOffTrials2 = OptoOffTrialsTemp(~isnan(OptoOffTrialsTemp));
OptoOnTrials = sData.behavior.opto.OptoOnTrialsIndices;
OptoAfterTrials = sData.behavior.opto.AfterOptoTrialsIndices;
%}

sData.imdata.gainMod.InOutRatio = struct;

nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.PlaceCellinAnyProt;
PFStartBin = NaN(nROIs,1);
PFLength = NaN(nROIs,1);
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
                
% generate a common PlaceFieldStartBin and Length Array for the three protocols (some cells are place cells in one, but not in another protocol. I want to include them)
for i = 1:1:nROIs
    if sum(ismember(PlaceCells,i))== 0
        continue
    end
    % Which place field is the largest (which bin, which protocol) , search for peak in position tuning curve
    % Opto-off
    PlaceFieldStartBinOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOff = PlaceFieldStartBinOffPre(~isnan(PlaceFieldStartBinOffPre));
    PlaceFieldLengthOffPre = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOff = PlaceFieldLengthOffPre(~isnan(PlaceFieldLengthOffPre));
    Peak = NaN(length(PlaceFieldLengthOff),1);
    for m = 1:1:length(PlaceFieldLengthOff) % if more then on place fields occurs
        Peak(m) = max(PosTuningOff(i,PlaceFieldStartBinOff(m):PlaceFieldStartBinOff(m)+PlaceFieldLengthOff(m)-1));
    end
    PeakOff = max(Peak);
    PeakPosOff = PlaceFieldStartBinOff(find(PeakOff));
    % Opto-on
    PlaceFieldStartBinOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOn = PlaceFieldStartBinOnPre(~isnan(PlaceFieldStartBinOnPre));
    PlaceFieldLengthOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOn = PlaceFieldLengthOnPre(~isnan(PlaceFieldLengthOnPre));
    Peak = NaN(length(PlaceFieldStartBinOn),1);
    for m = 1:1:length(PlaceFieldStartBinOn) % if more then on place fields occurs
        Peak(m) = max(PosTuningOn(i,PlaceFieldStartBinOn(m):PlaceFieldStartBinOn(m)+PlaceFieldLengthOn(m)-1));
    end
    PeakOn = max(Peak);
    PeakPosOn = PlaceFieldStartBinOn(find(PeakOn));
    % Opto-after
    PlaceFieldStartBinAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinAfter = PlaceFieldStartBinAfterPre(~isnan(PlaceFieldStartBinAfterPre));
    PlaceFieldLengthAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,:);
    PlaceFieldLengthAfter = PlaceFieldLengthAfterPre(~isnan(PlaceFieldLengthAfterPre));
    Peak = NaN(length(PlaceFieldStartBinAfter),1);
    for m = 1:1:length(PlaceFieldStartBinAfter) % if more then on place fields occurs
        Peak(m) = max(PosTuningAfter(i,PlaceFieldStartBinAfter(m):PlaceFieldStartBinAfter(m)+PlaceFieldLengthAfter(m)-1));
    end
    PeakAfter = max(Peak);
    PeakPosAfter = PlaceFieldStartBinAfter(find(PeakAfter));
    % compare place field peaks
    OffOnAfterPeak = [PeakOn PeakOff PeakAfter];
    if max(OffOnAfterPeak) == PeakOff
        PFStartBin(i) = PeakPosOff;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,find(PlaceFieldStartBinOffPre==PFStartBin(i))); 
    elseif max(OffOnAfterPeak) == PeakOn
        PFStartBin(i) = PeakPosOn;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,find(PlaceFieldStartBinOnPre==PFStartBin(i))); 
    elseif max(OffOnAfterPeak) == PeakAfter    
        PFStartBin(i) = PeakPosAfter;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,find(PlaceFieldStartBinAfterPre==PFStartBin(i))); 
    end
end
sData.imdata.gainMod.PFStartBin = PFStartBin;
sData.imdata.gainMod.PFLength = PFLength;
sData.imdata.gainMod.note = 'PFStartBin and Length: Place cells in any opto protocol, searched for largest peak among PFs';

% InOutRatio
for k = 1:1:3
     ActInPlaceField = NaN(nROIs,1);
     ActOutPlaceField = NaN(nROIs,1);
     ActInOutRatio = NaN(nROIs,1);
     if k == 1
         PosTuning = PosTuningOff;
     elseif k == 2
         PosTuning = PosTuningOn;
     elseif k == 3 
         PosTuning = PosTuningAfter;
     end
    % repeat for all three protocols
    % calculate mean-in-field/mean-out-of-field activity in opto-off vs opto-on trials
    for j = 1:1:length(PlaceCells)
        i = PlaceCells(j);
        ActInPlaceField(i) = mean(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));  
        ActOutPlaceField(i) = mean(PosTuning(i,PFStartBin(i)+PFLength(i):(PFStartBin(i)+nBins)-1));  
        ActInOutRatio(i) = ActInPlaceField(i)/ActOutPlaceField(i);
    end
    if k == 1
        sData.imdata.gainMod.InOutRatio.OptoOff.ActInPlaceField = ActInPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoOff.ActOutPlaceField = ActOutPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoOff.ActInOutRatio = ActInOutRatio;
    elseif k == 2
        sData.imdata.gainMod.InOutRatio.OptoOn.ActInPlaceField = ActInPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoOn.ActOutPlaceField = ActOutPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoOn.ActInOutRatio = ActInOutRatio;
    elseif k == 3
        sData.imdata.gainMod.InOutRatio.OptoAfter.ActInPlaceField = ActInPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoAfter.ActOutPlaceField = ActOutPlaceField;
        sData.imdata.gainMod.InOutRatio.OptoAfter.ActInOutRatio = ActInOutRatio;
    end
end

% generate summary
sData.imdata.gainMod.InOutRatio.ActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.imdata.gainMod.InOutRatio.ActInPlaceField_OffOnAfter(1:nROIs,1) = sData.imdata.gainMod.InOutRatio.OptoOff.ActInPlaceField;
sData.imdata.gainMod.InOutRatio.ActInPlaceField_OffOnAfter(1:nROIs,2) = sData.imdata.gainMod.InOutRatio.OptoOn.ActInPlaceField;
sData.imdata.gainMod.InOutRatio.ActInPlaceField_OffOnAfter(1:nROIs,3) = sData.imdata.gainMod.InOutRatio.OptoAfter.ActInPlaceField;

sData.imdata.gainMod.InOutRatio.ActOutPlaceField_OffOnAfter = NaN(nROIs,3);
sData.imdata.gainMod.InOutRatio.ActOutPlaceField_OffOnAfter(1:nROIs,1) = sData.imdata.gainMod.InOutRatio.OptoOff.ActOutPlaceField;
sData.imdata.gainMod.InOutRatio.ActOutPlaceField_OffOnAfter(1:nROIs,2) = sData.imdata.gainMod.InOutRatio.OptoOn.ActOutPlaceField;
sData.imdata.gainMod.InOutRatio.ActOutPlaceField_OffOnAfter(1:nROIs,3) = sData.imdata.gainMod.InOutRatio.OptoAfter.ActOutPlaceField;


sData.imdata.gainMod.InOutRatio.ActInOutRatio_OffOnAfter = NaN(nROIs,3);
sData.imdata.gainMod.InOutRatio.ActInOutRatio_OffOnAfter(1:nROIs,1) = sData.imdata.gainMod.InOutRatio.OptoOff.ActInOutRatio;
sData.imdata.gainMod.InOutRatio.ActInOutRatio_OffOnAfter(1:nROIs,2) = sData.imdata.gainMod.InOutRatio.OptoOn.ActInOutRatio;
sData.imdata.gainMod.InOutRatio.ActInOutRatio_OffOnAfter(1:nROIs,3) = sData.imdata.gainMod.InOutRatio.OptoAfter.ActInOutRatio;

sData.imdata.gainMod.InOutRatio.meanActInOutRatio_OptoOff = nanmean(sData.imdata.gainMod.InOutRatio.OptoOff.ActInOutRatio);
sData.imdata.gainMod.InOutRatio.meanActInOutRatio_OptoOn = nanmean(sData.imdata.gainMod.InOutRatio.OptoOn.ActInOutRatio);
sData.imdata.gainMod.InOutRatio.meanActInOutRatio_OptoAfter = nanmean(sData.imdata.gainMod.InOutRatio.OptoAfter.ActInOutRatio);

sData.imdata.gainMod.InOutRatio.stdActInOutRatio_OptoOff = nanstd(sData.imdata.gainMod.InOutRatio.OptoOff.ActInOutRatio);
sData.imdata.gainMod.InOutRatio.stdActInOutRatio_OptoOn = nanstd(sData.imdata.gainMod.InOutRatio.OptoOn.ActInOutRatio);
sData.imdata.gainMod.InOutRatio.stdActInOutRatio_OptoAfter = nanstd(sData.imdata.gainMod.InOutRatio.OptoAfter.ActInOutRatio);


end

