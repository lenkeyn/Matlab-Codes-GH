function sData = EffectPCinCtr(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects all place cells which are categorized as a place cell in control.
% compare peak amlitude and position of peak in control and opto trials
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols
% option to discard data before bin 10, also discard those place cells
% which has peak in the first 10 bin or categorized as landmark cells

%%% CALCULATE EFFECT SIZE OF THE OPTO STIM ON THE MEAN POS TUNING CURVES

sData.EffectPCinCtr = struct;

BinsToDiscard = 10; % discard data in the beginning of track
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells;
PlaceCellsPeakAfterBin10Pre = PlaceCells;
LandmarkCells = sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.LandmarkCellsList;

PFStartBin = NaN(nROIs,1);
PFLength = NaN(nROIs,1);
PFPeakAmpl = NaN(nROIs,1);
PFPeakPos = NaN(nROIs,1);
PosTuningOff = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig]; % generate circular data
PosTuningOn = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig]; % generate circular data
PosTuningAfter = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig]; % generate circular data
                
% generate a common PFStartBin and Length Array for the three protocols (some cells are place cells in one, but not in another protocol. I want to include them)
for j = 1:1:length(PlaceCells)
    i = PlaceCells(j); % roi ID of place cells
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
    PeakStartOff = PlaceFieldStartBinOff(find(PeakOff));
    PeakPosOff = min(find(PosTuningOff(i,1:sData.behavior.meta.nBins)==PeakOff));
    % Opto-on
    PlaceFieldStartBinOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinOn = PlaceFieldStartBinOnPre(~isnan(PlaceFieldStartBinOnPre));
    PlaceFieldLengthOnPre = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,:);
    PlaceFieldLengthOn = PlaceFieldLengthOnPre(~isnan(PlaceFieldLengthOnPre));
    if sum(~isnan(PlaceFieldStartBinOnPre))==0 % If there is no peak in the pos tuning curve
        PeakOn = 0;
    else        
        Peak = NaN(length(PlaceFieldStartBinOn),1);
        for m = 1:1:length(PlaceFieldStartBinOn) % if more then on place fields occurs
            Peak(m) = max(PosTuningOn(i,PlaceFieldStartBinOn(m):PlaceFieldStartBinOn(m)+PlaceFieldLengthOn(m)-1));
        end
        PeakOn = max(Peak);
        PeakStartOn = PlaceFieldStartBinOn(find(PeakOn));
        PeakPosOn = min(find(PosTuningOn(i,1:sData.behavior.meta.nBins)==PeakOn)); 
    end
    % Opto-after
    PlaceFieldStartBinAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldStartBin(i,:);
    PlaceFieldStartBinAfter = PlaceFieldStartBinAfterPre(~isnan(PlaceFieldStartBinAfterPre));
    PlaceFieldLengthAfterPre = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,:);
    PlaceFieldLengthAfter = PlaceFieldLengthAfterPre(~isnan(PlaceFieldLengthAfterPre));
    if sum(~isnan(PlaceFieldStartBinAfterPre))==0 % If there is no peak in the pos tuning curve
        PeakAfter = 0;
    else
        Peak = NaN(length(PlaceFieldStartBinAfter),1);
        for m = 1:1:length(PlaceFieldStartBinAfter) % if more then on place fields occurs
            Peak(m) = max(PosTuningAfter(i,PlaceFieldStartBinAfter(m):PlaceFieldStartBinAfter(m)+PlaceFieldLengthAfter(m)-1));
        end
        PeakAfter = max(Peak);
        PeakStartAfter = PlaceFieldStartBinAfter(find(PeakAfter));
        PeakPosAfter = min(find(PosTuningAfter(i,1:sData.behavior.meta.nBins)==PeakAfter)); 
    end
    % compare place field peaks
    OffOnAfterPeak = [PeakOn PeakOff PeakAfter];
    if max(OffOnAfterPeak) == PeakOff
        PFStartBin(i) = PeakStartOff;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,find(PlaceFieldStartBinOffPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakOff;
        PFPeakPos(i) = PeakPosOff;
    elseif max(OffOnAfterPeak) == PeakOn
        PFStartBin(i) = PeakStartOn;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceFieldBinLength(i,find(PlaceFieldStartBinOnPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakOn;
        PFPeakPos(i) = PeakPosOn;
    elseif max(OffOnAfterPeak) == PeakAfter    
        PFStartBin(i) = PeakStartAfter;
        PFLength(i) = sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceFieldBinLength(i,find(PlaceFieldStartBinAfterPre==PFStartBin(i))); 
        PFPeakAmpl(i) = PeakAfter;
        PFPeakPos(i) = PeakPosAfter;
    end
    if PFPeakPos(i) <= BinsToDiscard
        PlaceCellsPeakAfterBin10Pre(j) = NaN;
    end
end
for i = 1:1:length(LandmarkCells)
    PlaceCellsPeakAfterBin10Pre(PlaceCellsPeakAfterBin10Pre==(LandmarkCells(i)))= NaN;
end
PlaceCellsPeakAfterBin10 = PlaceCellsPeakAfterBin10Pre(~isnan(PlaceCellsPeakAfterBin10Pre));
sData.EffectPCinCtr.PlaceCellsPeakAfterBin10 = PlaceCellsPeakAfterBin10;
sData.EffectPCinCtr.PFStartBin = PFStartBin;
sData.EffectPCinCtr.PFLength = PFLength;
sData.EffectPCinCtr.PFPeakAmpl = PFPeakAmpl;
sData.EffectPCinCtr.PFPeakPos = PFPeakPos;
sData.EffectPCinCtr.note = 'PFStartBin Length Peak Ampl and Pos: searched for largest peak among PFs in all protocols, use the largest peak as place field peak and pos';

% Calculations
for k = 1:1:3
     meanActInPlaceField = NaN(nROIs,1); % mean Ca activity in the PF  
     peakActInPlaceField = NaN(nROIs,1); % peak amplitude in PF (mean)
     meanActOutPlaceField = NaN(nROIs,1); % mean Ca activity outside of the PF
     meanActInOutActRatio = NaN(nROIs,1); % ratio of mean activity inside/outside the PF
     peakActInMeanOutRatio = NaN(nROIs,1); % ratio of peak activity inside/ mean activity outside the PF
     BinPosPeakActInPlaceField = NaN(nROIs,1); % position (bin) of peak activity within PF (mean)
     %
     nROIsSelected = length(PlaceCellsPeakAfterBin10);
     PC2meanActInPlaceField = NaN(nROIsSelected,1); % mean Ca activity in the PF  
     PC2peakActInPlaceField = NaN(nROIsSelected,1); % peak amplitude in PF (mean)
     PC2meanActOutPlaceField = NaN(nROIsSelected,1); % mean Ca activity outside of the PF
     PC2meanActInOutActRatio = NaN(nROIsSelected,1); % ratio of mean activity inside/outside the PF
     PC2peakActInMeanOutRatio = NaN(nROIsSelected,1); % ratio of peak activity inside/ mean activity outside the PF
     PC2BinPosPeakActInPlaceField = NaN(nROIsSelected,1);
     if k == 1
         PosTuning = PosTuningOff;
         PosTuningWo10Bins = PosTuningOff;
         PosTuningWo10Bins(1:BinsToDiscard) = NaN;
         PosTuningWo10Bins(1:sData.behavior.meta.nBins + BinsToDiscard) = NaN;
     elseif k == 2
         PosTuning = PosTuningOn;
         PosTuningWo10Bins = PosTuningOn;
         PosTuningWo10Bins(1:BinsToDiscard) = NaN;
         PosTuningWo10Bins(1:sData.behavior.meta.nBins + BinsToDiscard) = NaN;
     elseif k == 3 
         PosTuning = PosTuningAfter;
         PosTuningWo10Bins = PosTuningAfter;
         PosTuningWo10Bins(1:BinsToDiscard) = NaN;
         PosTuningWo10Bins(1:sData.behavior.meta.nBins + BinsToDiscard) = NaN;
     end
    % repeat for all three protocols
    % calculate mean-in-field/mean-out-of-field activity in opto-off vs opto-on trials
    for j = 1:1:length(PlaceCells)
        i = PlaceCells(j);
        meanActInPlaceField(i) = mean(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));  
        peakActInPlaceField(i) = max(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));
        BinPosPeakActInPlaceField(i) = find(PosTuning(i,:)==peakActInPlaceField(i),1);
        meanActOutPlaceField(i) = mean(PosTuning(i,PFStartBin(i)+PFLength(i):(PFStartBin(i)+nBins)-1));  
        meanActInOutActRatio(i) = meanActInPlaceField(i)/meanActOutPlaceField(i);
        peakActInMeanOutRatio(i) = peakActInPlaceField(i)/meanActOutPlaceField(i);
    end
    for j = 1:1:length(PlaceCellsPeakAfterBin10) %selection of place cells: no landmark cells, cells discarded having peak before bin 10
        i = PlaceCellsPeakAfterBin10(j);
        PC2meanActInPlaceField(i) = mean(PosTuningWo10Bins(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));  
        PC2peakActInPlaceField(i) = max(PosTuningWo10Bins(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));
        PC2BinPosPeakActInPlaceField(i) = find(PosTuningWo10Bins(i,:)==peakActInPlaceField(i),1);
        PC2meanActOutPlaceField(i) = nanmean(PosTuningWo10Bins(i,PFStartBin(i)+PFLength(i):(PFStartBin(i)+nBins)-1));  
        PC2meanActInOutActRatio(i) = PC2meanActInPlaceField(i)/PC2meanActOutPlaceField(i);
        PC2peakActInMeanOutRatio(i) = PC2peakActInPlaceField(i)/PC2meanActOutPlaceField(i);
    end
    % organize data
    if k == 1
        sData.EffectPCinCtr.AllPC.OptoOff.meanActInPF = meanActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOff.peakActInPF = peakActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOff.BinPosPeakActInPF = BinPosPeakActInPlaceField; 
        sData.EffectPCinCtr.AllPC.OptoOff.meanActOutPF = meanActOutPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOff.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.EffectPCinCtr.AllPC.OptoOff.peakInMeanOutRatio = peakActInMeanOutRatio;
        % place cell selection
        sData.EffectPCinCtr.PC2.OptoOff.meanActInPF = PC2meanActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoOff.peakActInPF = PC2peakActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoOff.BinPosPeakActInPF = PC2BinPosPeakActInPlaceField; 
        sData.EffectPCinCtr.PC2.OptoOff.meanActOutPF = PC2meanActOutPlaceField;
        sData.EffectPCinCtr.PC2.OptoOff.meanActInOutMeanActRatio = PC2meanActInOutActRatio;
        sData.EffectPCinCtr.PC2.OptoOff.peakInMeanOutRatio = PC2peakActInMeanOutRatio;
    elseif k == 2
        sData.EffectPCinCtr.AllPC.OptoOn.meanActInPF = meanActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOn.peakActInPF = peakActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOn.BinPosPeakActInPF = BinPosPeakActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOn.meanActOutPF = meanActOutPlaceField;
        sData.EffectPCinCtr.AllPC.OptoOn.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.EffectPCinCtr.AllPC.OptoOn.peakInMeanOutRatio = peakActInMeanOutRatio;
        % selected place cells
        sData.EffectPCinCtr.PC2.OptoOn.meanActInPF = PC2meanActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoOn.peakActInPF = PC2peakActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoOn.BinPosPeakActInPF = PC2BinPosPeakActInPlaceField; 
        sData.EffectPCinCtr.PC2.OptoOn.meanActOutPF = PC2meanActOutPlaceField;
        sData.EffectPCinCtr.PC2.OptoOn.meanActInOutMeanActRatio = PC2meanActInOutActRatio;
        sData.EffectPCinCtr.PC2.OptoOn.peakInMeanOutRatio = PC2peakActInMeanOutRatio;
    elseif k == 3
        sData.EffectPCinCtr.AllPC.OptoAfter.meanActInPF = meanActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoAfter.peakActInPF = peakActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoAfter.BinPosPeakActInPF = BinPosPeakActInPlaceField;
        sData.EffectPCinCtr.AllPC.OptoAfter.meanActOutPF = meanActOutPlaceField;
        sData.EffectPCinCtr.AllPC.OptoAfter.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.EffectPCinCtr.AllPC.OptoAfter.peakInMeanOutRatio = peakActInMeanOutRatio;
        % selected place cells
        sData.EffectPCinCtr.PC2.OptoAfter.meanActInPF = PC2meanActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoAfter.peakActInPF = PC2peakActInPlaceField;
        sData.EffectPCinCtr.PC2.OptoAfter.BinPosPeakActInPF = PC2BinPosPeakActInPlaceField; 
        sData.EffectPCinCtr.PC2.OptoAfter.meanActOutPF = PC2meanActOutPlaceField;
        sData.EffectPCinCtr.PC2.OptoAfter.meanActInOutMeanActRatio = PC2meanActInOutActRatio;
        sData.EffectPCinCtr.PC2.OptoAfter.peakInMeanOutRatio = PC2peakActInMeanOutRatio;
    end
end

% generate summary ALL PCs
nPC = length(PlaceCells);
sData.EffectPCinCtr.AllPC.AllPCandLandmarkCells = PlaceCells;
sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter = NaN(nPC,3);
sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.meanActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.meanActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.meanActInPF(PlaceCells);

sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter = NaN(nPC,3);
sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.peakActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.peakActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.peakActInPF(PlaceCells);

sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter = NaN(nPC,3);
sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.BinPosPeakActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.BinPosPeakActInPF(PlaceCells);
sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.BinPosPeakActInPF(PlaceCells);

sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter = NaN(nPC,3);
sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.meanActOutPF(PlaceCells);
sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.meanActOutPF(PlaceCells);
sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.meanActOutPF(PlaceCells);

sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter = NaN(nPC,3);
sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.meanActInOutMeanActRatio(PlaceCells);
sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.meanActInOutMeanActRatio(PlaceCells);
sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.meanActInOutMeanActRatio(PlaceCells);

sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter = NaN(nPC,3); 
sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(1:nPC,1) = sData.EffectPCinCtr.AllPC.OptoOff.peakInMeanOutRatio(PlaceCells);
sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(1:nPC,2) = sData.EffectPCinCtr.AllPC.OptoOn.peakInMeanOutRatio(PlaceCells);
sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(1:nPC,3) = sData.EffectPCinCtr.AllPC.OptoAfter.peakInMeanOutRatio(PlaceCells);


% means
sData.EffectPCinCtr.AllPC.means = struct;
sData.EffectPCinCtr.AllPC.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.AllPC.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.AllPC.meanInMeanOutRatio_OffOnAfter(:,3));

sData.EffectPCinCtr.AllPC.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.AllPC.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.AllPC.peakInMeanOutRatio_OffOnAfter(:,3));

sData.EffectPCinCtr.AllPC.means.peakActInPF_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.AllPC.means.peakActInPF_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.peakActInPF_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.peakActInPF_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.AllPC.peakActInPF_OffOnAfter(:,3));

sData.EffectPCinCtr.AllPC.means.BinPosPeakActInPF_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.AllPC.means.BinPosPeakActInPF_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.BinPosPeakActInPF_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.BinPosPeakActInPF_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.AllPC.BinPosPeakActInPF_OffOnAfter(:,2));

sData.EffectPCinCtr.AllPC.means.InOff = nanmean(sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.InOn= nanmean(sData.EffectPCinCtr.AllPC.meanActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.OutOff = nanmean(sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter(:,1));
sData.EffectPCinCtr.AllPC.means.OutOn = nanmean(sData.EffectPCinCtr.AllPC.meanActOutPF_OffOnAfter(:,2));
sData.EffectPCinCtr.AllPC.means.InOnPerInOff = sData.EffectPCinCtr.AllPC.means.InOn/sData.EffectPCinCtr.AllPC.means.InOff;
sData.EffectPCinCtr.AllPC.means.OutOnPerOutOff = sData.EffectPCinCtr.AllPC.means.OutOn/sData.EffectPCinCtr.AllPC.means.OutOff;

% generate summary selected PCs
nROIsSelected = length(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.PlaceCellSelection = PlaceCellsPeakAfterBin10;
sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter = NaN(nROIsSelected,3);
sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.meanActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.meanActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.meanActInPF(PlaceCellsPeakAfterBin10);

sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter = NaN(nROIsSelected,3);
sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.peakActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.peakActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.peakActInPF(PlaceCellsPeakAfterBin10);

sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter = NaN(nROIsSelected,3);
sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.BinPosPeakActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.BinPosPeakActInPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.BinPosPeakActInPF(PlaceCellsPeakAfterBin10);

sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter = NaN(nROIsSelected,3);
sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.meanActOutPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.meanActOutPF(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.meanActOutPF(PlaceCellsPeakAfterBin10);

sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter = NaN(nROIsSelected,3);
sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.meanActInOutMeanActRatio(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.meanActInOutMeanActRatio(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.meanActInOutMeanActRatio(PlaceCellsPeakAfterBin10);

sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter = NaN(nROIsSelected,3); 
sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(1:nROIsSelected,1) = sData.EffectPCinCtr.PC2.OptoOff.peakInMeanOutRatio(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(1:nROIsSelected,2) = sData.EffectPCinCtr.PC2.OptoOn.peakInMeanOutRatio(PlaceCellsPeakAfterBin10);
sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(1:nROIsSelected,3) = sData.EffectPCinCtr.PC2.OptoAfter.peakInMeanOutRatio(PlaceCellsPeakAfterBin10);

% means
sData.EffectPCinCtr.PC2.means = struct;
sData.EffectPCinCtr.PC2.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.PC2.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.PC2.meanInMeanOutRatio_OffOnAfter(:,3));

sData.EffectPCinCtr.PC2.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.PC2.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.PC2.peakInMeanOutRatio_OffOnAfter(:,3));

sData.EffectPCinCtr.PC2.means.peakActInPF_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.PC2.means.peakActInPF_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.peakActInPF_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.peakActInPF_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,3));

sData.EffectPCinCtr.PC2.means.BinPosPeakActInPF_OffOnAfter = NaN(1,3);
sData.EffectPCinCtr.PC2.means.BinPosPeakActInPF_OffOnAfter(1,1) = nanmean(sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.BinPosPeakActInPF_OffOnAfter(1,2) = nanmean(sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.BinPosPeakActInPF_OffOnAfter(1,3) = nanmean(sData.EffectPCinCtr.PC2.BinPosPeakActInPF_OffOnAfter(:,2));

sData.EffectPCinCtr.PC2.means.InOff = nanmean(sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.InOn= nanmean(sData.EffectPCinCtr.PC2.meanActInPF_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.OutOff = nanmean(sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter(:,1));
sData.EffectPCinCtr.PC2.means.OutOn = nanmean(sData.EffectPCinCtr.PC2.meanActOutPF_OffOnAfter(:,2));
sData.EffectPCinCtr.PC2.means.InOnPerInOff = sData.EffectPCinCtr.PC2.means.InOn/sData.EffectPCinCtr.PC2.means.InOff;
sData.EffectPCinCtr.PC2.means.OutOnPerOutOff = sData.EffectPCinCtr.PC2.means.OutOn/sData.EffectPCinCtr.PC2.means.OutOff;
sData.EffectPCinCtr.PC2.note = 'All place cells which are not landmark cell and their peak is not in the first 10 bins';

%%% CALCULATE IF THE EFFECT WAS SIGNIFICANT ON THE CELL LEVEL
% use non parametric rank sum test on the peak in the PF in ctr and opto
sData.EffectPCinCtr.ranksumCellEffect = struct;
ranksumP = NaN(nROIsSelected,1); % P value for the Unsigned Wilcokon test
ranksumH = NaN(nROIsSelected,1); % H=0 non-sign, H=1 sign indicates that the null hypothesis can be rejected 
OptoOnTrials = sData.behavior.opto.OptoOnTrialsIndices;
OptoOffTrials = sData.behavior.opto.OptoOffTrialsIndices;
%PeaksInPFOff3bin = NaN(length(sData.behavior.opto.OptoOffTrialsIndices),3);
%PeaksInPFOn3bin = NaN(length(sData.behavior.opto.OptoOnTrialsIndices),3);
CellEffectPeak = NaN(nROIsSelected,1);
nonsignificant = 0;
signIncrease = 0;
signDecrease = 0;
for i = 1:1:nROIsSelected
    j = PlaceCellsPeakAfterBin10(i);
    PFPeakPos = sData.EffectPCinCtr.PFPeakPos(j);
    data2 = horzcat(sData.imdata.binned.RoidFF{1, j},sData.imdata.binned.RoidFF{1, j}); %binned Ca activity for each roi
    PeaksInPFOff3bin = data2(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1);
    PeaksInPFOff = mean(PeaksInPFOff3bin,2);
    PeaksInPFOn3bin = data2(OptoOnTrials,PFPeakPos-1:PFPeakPos+1); % data(OptoOffTrials,PFPeakPos-1:PFPeakPos+1); use 1 or 3 bins?
    PeaksInPFOn = mean(PeaksInPFOn3bin,2);
    [ranksumP(i),ranksumH(i)] = ranksum(PeaksInPFOff,PeaksInPFOn); % h=1 significant, h=0 not sign
    CellEffectPeak(i) = mean(PeaksInPFOn)/mean(PeaksInPFOff); % less than 1 decrease
    if ranksumH(i) == 0
       nonsignificant = nonsignificant + 1; 
    elseif CellEffectPeak(i) > 1 && ranksumH(i) == 1
       signIncrease = signIncrease + 1;
    elseif CellEffectPeak(i) < 1 && ranksumH(i) == 1
       signDecrease = signDecrease + 1;
    end
end

sData.EffectPCinCtr.ranksumCellEffect.ROIs = PlaceCellsPeakAfterBin10;
sData.EffectPCinCtr.ranksumCellEffect.CellEffectPeakOnPerOff = CellEffectPeak;
sData.EffectPCinCtr.ranksumCellEffect.ranksumP = ranksumP;
sData.EffectPCinCtr.ranksumCellEffect.ranksumH = ranksumH;
sData.EffectPCinCtr.ranksumCellEffect.SignIncrNu = signIncrease;
sData.EffectPCinCtr.ranksumCellEffect.SignDecrNu = signDecrease;
sData.EffectPCinCtr.ranksumCellEffect.nonsignificantNu = nonsignificant;
sData.EffectPCinCtr.ranksumCellEffect.SignIncrPercentage = signIncrease/nROIsSelected*100;
sData.EffectPCinCtr.ranksumCellEffect.SignDecrPercentage = signDecrease/nROIsSelected*100;
sData.EffectPCinCtr.ranksumCellEffect.nonsignificantPercentage = nonsignificant/nROIsSelected*100;


% added 2023.06.30.
%%% calculate percentage of increase and decrease on cell level (no matter if it was significant or not)
Incr = 0; Decr = 0;
PeakOff = sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,1);
PeakOn = sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,2);
EffectRatio = PeakOff.\PeakOn;
for i = 1:1:nROIsSelected
    if EffectRatio(i) > 1 
        Incr = Incr + 1;
    elseif EffectRatio(i) < 1
        Decr = Decr + 1;
    end
end
PercentageCellEffectIncr = Incr / nROIsSelected;
PercentageCellEffectDecr = Decr / nROIsSelected;
sData.EffectPCinCtr.PercCellEffectIncr = PercentageCellEffectIncr;
sData.EffectPCinCtr.PercCellEffectDecr = PercentageCellEffectDecr;

figure('Color','white')
pie([PercentageCellEffectDecr PercentageCellEffectIncr]);


%%% CALCULATE IF THE EFFECT WAS SIGNIFICANT ON SESSION LEVEL % added 2023.06.30.
% use non parametric PAIRED rank sum test (Wilcoxon) on the peak in the PF
% in ctr and opto for each cell
PeakOff = sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,1);
PeakOn = sData.EffectPCinCtr.PC2.peakActInPF_OffOnAfter(:,2);
% WPvalue: P value for the Unsigned Wilcokon test
% WH: H=0 non-sign, H=1 sign indicates that the null hypothesis can be rejected 

%CellEffectPeak = NaN(nROIsSelected,1);
[WPvalue,WH] = signrank(PeakOff,PeakOn); 
CellEffectPeak = mean(PeakOn)/mean(PeakOff); % less than 1 decrease

sData.EffectPCinCtr.WilcoxonEffectSession.CellEffectPeakInPF = CellEffectPeak;
sData.EffectPCinCtr.WilcoxonEffectSession.Pvalue = WPvalue;
sData.EffectPCinCtr.WilcoxonEffectSession.ifSignificant = WH;


%%% CALCULATE ZERO SHIFTED POS TUNING CURVES AND CALCULATE SHIFT AND GAIN FOR EACH CELL
sData.EffectPCinCtr.PosTuningShiftedCells = struct;
BinsAroundPeak = 30;
PosTuningShiftedOff = NaN(nROIsSelected,2*BinsAroundPeak+1);
PosTuningShiftedOn = NaN(nROIsSelected,2*BinsAroundPeak+1);
ShiftEffectSize = NaN(nROIsSelected,1);
GainEffectSize = NaN(nROIsSelected,1);
GainShiftRatio = NaN(nROIsSelected,1);
GainEffectSizePercentage = NaN(nROIsSelected,1);
ShiftEffectSizePercentage = NaN(nROIsSelected,1);

for i = 1:1:nROIsSelected
    j = PlaceCellsPeakAfterBin10(i);
    PFPeakPos = sData.EffectPCinCtr.PFPeakPos(j);
    PTOff = PosTuningOff(j,:);
    PTOn = PosTuningOn(j,:);
    if PFPeakPos > BinsAroundPeak
        PosTuningShiftedOff(i,:) = PTOff((PFPeakPos-BinsAroundPeak):(PFPeakPos+BinsAroundPeak));
        PosTuningShiftedOn(i,:) = PTOn((PFPeakPos-BinsAroundPeak):(PFPeakPos+BinsAroundPeak));
    else
        PosTuningShiftedOff(i,:) = PTOff((PFPeakPos+sData.behavior.meta.nBins-BinsAroundPeak):(PFPeakPos+sData.behavior.meta.nBins+BinsAroundPeak));
        PosTuningShiftedOn(i,:) = PTOn((PFPeakPos+sData.behavior.meta.nBins-BinsAroundPeak):(PFPeakPos+sData.behavior.meta.nBins+BinsAroundPeak));
    end
    % gain moulation and shift calculation for each Place cell ROIs
    AmplChangeinPFPeakOnMinusOff = PosTuningShiftedOn(i,BinsAroundPeak+1) - PosTuningShiftedOff(i,BinsAroundPeak+1);
    AmplChangeBeforePFOnMinusOff = PosTuningShiftedOn(i,BinsAroundPeak-24) - PosTuningShiftedOff(i,BinsAroundPeak-24); %baseline change before the peak
    ShiftEffectSize(i) = abs(AmplChangeBeforePFOnMinusOff);
    GainEffectSize(i) = abs(AmplChangeinPFPeakOnMinusOff - AmplChangeBeforePFOnMinusOff);
    GainShiftRatio(i) = GainEffectSize(i)/ShiftEffectSize(i);
    GainEffectSizePercentage(i) = GainEffectSize(i)/(GainEffectSize(i) + ShiftEffectSize(i));
    ShiftEffectSizePercentage(i) = ShiftEffectSize(i)/(GainEffectSize(i) + ShiftEffectSize(i));
end

sData.EffectPCinCtr.PosTuningShiftedCells.PosTuningShiftedOff = PosTuningShiftedOff;
sData.EffectPCinCtr.PosTuningShiftedCells.PosTuningShiftedOn = PosTuningShiftedOn;
sData.EffectPCinCtr.PosTuningShiftedCells.ShiftEffectSize = ShiftEffectSize;
sData.EffectPCinCtr.PosTuningShiftedCells.GainEffectSize = GainEffectSize;
sData.EffectPCinCtr.PosTuningShiftedCells.GainShiftRatio = GainShiftRatio;
sData.EffectPCinCtr.PosTuningShiftedCells.GainEffectSizePercentage = GainEffectSizePercentage;
sData.EffectPCinCtr.PosTuningShiftedCells.ShiftEffectSizePercentage = ShiftEffectSizePercentage;


% Gain modulation and shift on mean pos tuning curve
% calculate change in peak and in baseline, Peak is in bin 31, Baseline calculation is at bin 6, 25 bin before peak (50 cm)
sData.EffectPCinCtr.meanPosTuningShifted = struct;

meanPosTuningShiftedOff = nanmean(PosTuningShiftedOff,1);
meanPosTuningShiftedOn = nanmean(PosTuningShiftedOn,1); 
meanAmplChangeinPFPeakOnMinusOff = meanPosTuningShiftedOn(BinsAroundPeak+1) - meanPosTuningShiftedOff(BinsAroundPeak+1);
meanAmplChangeBeforePFOnMinusOff = meanPosTuningShiftedOn(BinsAroundPeak-24) - meanPosTuningShiftedOff(BinsAroundPeak-24); %baseline change before the peak
ShiftEffectSize = abs(meanAmplChangeBeforePFOnMinusOff);
GainEffectSize = abs(meanAmplChangeinPFPeakOnMinusOff - meanAmplChangeBeforePFOnMinusOff);
GainShiftRatio = GainEffectSize/ShiftEffectSize;
GainEffectSizePercentage = GainEffectSize/(GainEffectSize + ShiftEffectSize);
ShiftEffectSizePercentage = ShiftEffectSize/(GainEffectSize + ShiftEffectSize);

sData.EffectPCinCtr.meanPosTuningShifted.PosTuningShiftedOff = PosTuningShiftedOff;
sData.EffectPCinCtr.meanPosTuningShifted.PosTuningShiftedOn = PosTuningShiftedOn;
sData.EffectPCinCtr.meanPosTuningShifted.meanPosTuninShiftedOff = meanPosTuningShiftedOff;    
sData.EffectPCinCtr.meanPosTuningShifted.meanPosTuninShiftedOn = meanPosTuningShiftedOn;   
sData.EffectPCinCtr.meanPosTuningShifted.meanAmplChangeinPFPeakOnvsOff = meanAmplChangeinPFPeakOnMinusOff; 
sData.EffectPCinCtr.meanPosTuningShifted.meanAmplChangeBeforePFPeakOnvsOff = meanAmplChangeBeforePFOnMinusOff;  
sData.EffectPCinCtr.meanPosTuningShifted.ShiftEffectSize = ShiftEffectSize;
sData.EffectPCinCtr.meanPosTuningShifted.GainEffectSize = GainEffectSize;
sData.EffectPCinCtr.meanPosTuningShifted.GainShiftRatio = GainShiftRatio;
sData.EffectPCinCtr.meanPosTuningShifted.GainEffectSizePercentage = GainEffectSizePercentage;
sData.EffectPCinCtr.meanPosTuningShifted.ShiftEffectSizePercentage = ShiftEffectSizePercentage;

%%% calculate Regression lines (linear fit) to two datapoints: 1:(outofPCctr,outofPFopto) 2:(inPFctr,inPFopto)
sData.EffectPCinCtr.RegrForTwoPoints = struct;

EqSteepness = NaN(nROIsSelected,1); % steepness of equation
EqShift = NaN(nROIsSelected,1); % shift of equation
% caclulate 2 point: baseline (ctr,opto) and peak(ctr,opto)
PosTuningOffNaN = PosTuningOff; 
PosTuningOnNaN = PosTuningOn;
PosTuningOffNaN(:,1:BinsToDiscard) = NaN; % discard data from first bins
PosTuningOnNaN(:,1:BinsToDiscard) = NaN;
PosTuningOffNaN(:,sData.behavior.meta.nBins+1:sData.behavior.meta.nBins+BinsToDiscard) = NaN; 
PosTuningOnNaN(:,sData.behavior.meta.nBins+1:sData.behavior.meta.nBins+BinsToDiscard) = NaN;
%REGRESSION FOR 2 POINTS: not to overrepresent baseline points
for i = 1:1:nROIsSelected
    roi = PlaceCellsPeakAfterBin10(i);
    % activity at PF peak in ctr and opto
    PeakInPFOff3bin = nanmean(PosTuningOffNaN(roi,PFPeakPos-1:PFPeakPos+1)); 
    PeakInPFOn3bin = nanmean(PosTuningOnNaN(roi,PFPeakPos-1:PFPeakPos+1));
    % activity out of PF (baseline
    BaselineBins = NaN(6,1); % data of 6 bins will be averaged
    MinPosTuningOff = min(PosTuningOffNaN(roi,PFStartBin(roi)+PFLength(roi):sData.behavior.meta.nBins+PFStartBin(roi))); % search the minimum in the PosTuning curve expect the place fieald
    PosMinPosTuningOff = min(find(PosTuningOffNaN(roi,:) == MinPosTuningOff)); % position of minimum in opto-off
    BaselineBins(1:3) = [PosMinPosTuningOff-1 PosMinPosTuningOff PosMinPosTuningOff+1];
    MinPosTuningOn = min(PosTuningOnNaN(roi,PFStartBin(roi)+PFLength(roi):sData.behavior.meta.nBins+PFStartBin(roi)));
    PosMinPosTuningOn = min(find(PosTuningOnNaN(roi,:) == MinPosTuningOn)); % position of minimum in opto-on
    BaselineBins(4:6) = [PosMinPosTuningOn-1 PosMinPosTuningOn PosMinPosTuningOn+1];
    BBins = unique(BaselineBins);
    BaselineOutPFOff6Bins = nanmean(PosTuningOffNaN(roi,BBins));
    BaselineOutPFOn6Bins = nanmean(PosTuningOnNaN(roi,BBins));
    
    % calculate linear fit for data and plot the opto-off vs opto-on amplitude value for each roi for each bin
    PosTunOffTwoPoints = [BaselineOutPFOff6Bins BaselineOutPFOn6Bins];
    PosTunOnTwoPoints = [PeakInPFOff3bin PeakInPFOn3bin];
    %Max = max(max(PosTunOff),max(PosTunOnTwoPoints));
    [~, ~, R2, P, st, sh] = linFit(PosTunOffTwoPoints,PosTunOnTwoPoints);
    %Regression2(i) = R2;
    %SignificanceP(i) = P;
    EqSteepness(i) = st; % steepness of equation
    EqShift(i) = sh; % shift of equation
end
sData.EffectPCinCtr.RegrForTwoPoints.EqSteepness = EqSteepness; 
sData.EffectPCinCtr.RegrForTwoPoints.EqShift = EqShift; 
sData.EffectPCinCtr.RegrForTwoPoints.meanSteepness = mean(EqSteepness);
sData.EffectPCinCtr.RegrForTwoPoints.meanShift = mean(EqShift);
sData.EffectPCinCtr.RegrForTwoPoints.note = 'all place cell which has peak after bin 10 and not landmark cells, two points were generated for each cell: one ctr-opto pair where the PF peak (3 bin average) and 3+3 bin average of around the minimum in ctr and opto respectively ) is and '; 


%%% REGRESSION FOR ALL DATAPOINTS (PLACE CELL SELECTION
PosTuningOffOnArray = NaN(nROIsSelected*(sData.behavior.meta.nBins-BinsToDiscard),2);
PosTuningOffArrayPre = PosTuningOff(PlaceCellsPeakAfterBin10,BinsToDiscard+1:sData.behavior.meta.nBins);
PosTuningOnArrayPre = PosTuningOn(PlaceCellsPeakAfterBin10,BinsToDiscard+1:sData.behavior.meta.nBins);
PosTuningOffOnArray(:,1) = PosTuningOffArrayPre(:);
PosTuningOffOnArray(:,2) = PosTuningOnArrayPre(:);
[~, ~, R2, P, st, sh] = linFit(PosTuningOffOnArray(:,1),PosTuningOffOnArray(:,2));
sData.EffectPCinCtr.RegrAllDataPoints.PosTuningOffOnArray = PosTuningOffOnArray;
sData.EffectPCinCtr.RegrAllDataPoints.Regression2 = R2;
sData.EffectPCinCtr.RegrAllDataPoints.SignificanceP = P;
sData.EffectPCinCtr.RegrAllDataPoints.EqSteepness = st; 
sData.EffectPCinCtr.RegrAllDataPoints.EqShift = sh; 
sData.EffectPCinCtr.RegrAllDataPoints.note = 'all place cell which has peak after bin 10 and not landmark cells, all individual ctr-opto data pairs from the binned Ca data'; 

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end