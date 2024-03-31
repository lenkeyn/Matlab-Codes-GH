function sData = gainModulationInOutRatioPCinAllProtocols(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects all place cells which are categorized as a place cell by any protocol.
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols

sData.gainModulationPCinAllProtocol.InOutMeanActRatio = struct;

nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PlaceCells = sData.imdata.MaoPC_Opto_dff.PlaceCellinAllProt;
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
sData.gainModulationPCinAllProtocol.PFStartBin = PFStartBin;
sData.gainModulationPCinAllProtocol.PFLength = PFLength;
sData.gainModulationPCinAllProtocol.note = 'PFStartBin and Length: Place cells in any opto protocol, searched for largest peak among PFs';

% InOutMeanActRatio
for k = 1:1:3
     meanActInPlaceField = NaN(nROIs,1); % mean Ca activity in the PF  
     peakActInPlaceField = NaN(nROIs,1); % peak amplitude in PF (mean)
     meanActOutPlaceField = NaN(nROIs,1); % mean Ca activity outside of the PF
     meanActInOutActRatio = NaN(nROIs,1); % ratio of mean activity inside/outside the PF
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
        meanActInOutActRatio(i) = meanActInPlaceField(i)/meanActOutPlaceField(i);
        peakActInMeanOutRatio(i) = peakActInPlaceField(i)/meanActOutPlaceField(i);
    end
    if k == 1
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField; 
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 2
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 3
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio = peakActInMeanOutRatio;
    end
end

% generate summary
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActInPlaceField;

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.peakActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.peakActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.peakActInPlaceField;

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField;

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActOutPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActOutPlaceField;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActOutPlaceField;

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio;

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter = NaN(nROIs,3); 
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOff.peakInMeanOutRatio;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoOn.peakInMeanOutRatio;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAllProtocol.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio;


% means
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means = struct;
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,3));

sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinAllProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAllProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end