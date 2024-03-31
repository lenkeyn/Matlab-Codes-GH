function sData = gainModulationInOutRatioPCinAnyProtocols(sData)

% the function checks average postition tuning curves of ROIS without and with optical stimulation.
% collects Any place cells which are categorized as a place cell by any protocol.
% compare the mean Ca-activity within the place field and outside the place field
% compare peak Ca-activity in different opto protocols

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio = struct;

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
sData.gainModulationPCinAnyProtocol.PFStartBin = PFStartBin;
sData.gainModulationPCinAnyProtocol.PFLength = PFLength;
sData.gainModulationPCinAnyProtocol.note = 'PFStartBin and Length: Place cells in any opto protocol, searched for largest peak among PFs';

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
    % repeat for Any three protocols
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
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField; 
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 2
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.peakInMeanOutRatio = peakActInMeanOutRatio;
    elseif k == 3
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActInPlaceField = meanActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.peakActInPlaceField = peakActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField = BinPosPeakActInPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActOutPlaceField = meanActOutPlaceField;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio = meanActInOutActRatio;
        sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio = peakActInMeanOutRatio;
    end
end

% generate summary
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActInPlaceField;

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.peakActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.peakActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.peakActInPlaceField;

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.BinPosPeakActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.BinPosPeakActInPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.BinPosPeakActInPlaceField;

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActOutPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActOutPlaceField;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanActOutPlaceField_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActOutPlaceField;

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter = NaN(nROIs,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.meanActInOutMeanActRatio;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.meanActInOutMeanActRatio;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.meanActInOutMeanActRatio;

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter = NaN(nROIs,3); 
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,1) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOff.peakInMeanOutRatio;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,2) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoOn.peakInMeanOutRatio;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(1:nROIs,3) = sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.OptoAfter.peakInMeanOutRatio;


% means
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means = struct;
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.meanInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,1));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,2));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakInMeanOutRatio_OffOnAfter(:,3));

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.peakActInPlaceField_OffOnAfter(:,3));

sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter = NaN(1,3);
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,1));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));
sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3) = nanmean(sData.gainModulationPCinAnyProtocol.InOutMeanActRatio.BinPosPeakActInPlaceField_OffOnAfter(:,2));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end