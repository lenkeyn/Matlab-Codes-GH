function sData = gainModulationInOutRatioLandmarkAndOnePFPC(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects all place cells which are categorized as a place cell by any protocol.
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio = struct;

nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCellsPre = [sData.imdata.MaoPC_Opto_dff.LandmarkCells.LandmarkCellinAllProt; sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.PlaceCellsWithOnePFinAllProt];
PlaceCells = sort(PlaceCellsPre);
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
sData.gainModulationPCinLandmarkAndOnePFPC.PFStartBin = PFStartBin;
sData.gainModulationPCinLandmarkAndOnePFPC.PFLength = PFLength;
sData.gainModulationPCinLandmarkAndOnePFPC.note = 'PFStartBin and Length: landmark cells and place cells with only one place field in both protocols, searched for largest peak among PFs';

% InOutMeanActRatio
for k = 1:1:3
     meanActInPlaceField = NaN(nROIs,1); % mean Ca activity in the PF  
     peakActInPlaceField = NaN(nROIs,1); % peak amplitude in PF (mean)
     meanActOutPlaceField = NaN(nROIs,1); % mean Ca activity outside of the PF
     meanActInOutRatio = NaN(nROIs,1); % ratio of mean activity inside/outside the PF
     peakActInMeanOutRatio = NaN(nROIs,1); % ratio of peak activity inside/ mean activity outside the PF
     BinPosPeakActInPlaceField = NaN(nROIs,1); % position (bin) of peak activity within PF (mean)
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
        meanActInPlaceField(i) = mean(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));  
        peakActInPlaceField(i) = max(PosTuning(i,PFStartBin(i):(PFStartBin(i)+PFLength(i)-1)));
        BinPosPeakActInPlaceField(i) = find(PosTuning(i,:)==peakActInPlaceField(i),1);
        meanActOutPlaceField(i) = mean(PosTuning(i,PFStartBin(i)+PFLength(i):(PFStartBin(i)+nBins)-1));  
        meanActInOutRatio(i) = meanActInPlaceField(i)/meanActOutPlaceField(i);
        peakActInMeanOutRatio(i) = peakActInPlaceField(i)/meanActOutPlaceField(i);
    end
    if k == 1
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField; 
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio = meanActInOutRatio;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 2
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio = meanActInOutRatio;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 3
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio = meanActInOutRatio;
        sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio = peakActInMeanOutRatio;
    end
end

% generate summary
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActInPlaceField;

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.peakActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.peakActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.peakActInPlaceField;

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField;

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActOutPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActOutPlaceField;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActOutPlaceField;

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio;

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter = NaN(nROIs,3); 
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOff.peakInMeanOutRatio;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoOn.peakInMeanOutRatio;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio;


% means
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means = struct;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,3));

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));

sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOff = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOn= nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOff = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOn = nanmean(sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOnPerInOff = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOn/sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOff;
sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOnPerOutOff = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOn/sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOff;


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end