function sData = gainMod_RegressionForTwoPointsOpto(sData)

% REGRESSION CURVE WAS MADE FROM POSITION TUNING CURVE IN OPTO-ON vs OPTO-OFF conditions with script gainMod_binnedOptoOffOnRegression
% As I used 80 bins, I got 80 points for the regression curve (Opto-off vs % opto-on). 
% But the place field (PF) only takes usually 5-10 bins, so the not place-field positions dominate the regression. 
% Therefore I decided to derive one datapoint both from the PF-in and PF-out datapoints. 
% Selection of PF-in datapoint: search the maximum (bin) within place field (PF based on opto-off condtition).
% Selection of PF-out (baseline) datapoint: search for minimum Ca-signal in opto-off condition (except first bins of the trial, because those bins are discarded due to the delay of optical stimulation in some cases):

% create folder: GainModTwoPointRegression
% savePath = 'E:\ANALYSIS\m8061-VIPArch-Thy1-GCamp\m8061-20200523-00-second-anal\Imaging\GainModTwoPointRegression';
% mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'GainModTwoPointRegression');
% savePath = strcat(sData.sessionInfo.savePath,'\Imaging\GainModTwoPointRegression');

DiscardCmBeginning = 20; % discard the first bins, because optical stimulation sometimes delayed due to high speed of animal
% datasetDFF = sData.imdata.binned.RoidFF; % use dFF data
DiscardBins = round(DiscardCmBeginning/sData.behavior.meta.binSize);
PlaceCells = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells; % ROI-ID of place cells
nPC = length(PlaceCells); % number of place cells
nROIs = sData.imdata.nROIs;
nBins = sData.behavior.meta.nBins;
PeakPFin3BinsMean = NaN(nROIs,3); % columns: opto-off, opto-on, after-opto
BaselinePFout3BinsMean = NaN(nROIs,3); 
PeakPFin = NaN(nROIs,3);
PeakPFinPos = NaN(nROIs,1);
BaselinePFout = NaN(nROIs,3);
BaselinePFoutPos = NaN(nROIs,1);
BaselinePFout3percentile = NaN(nROIs,3); 
PFOutActScheme = NaN(nROIs,nBins);
LabelsForColumns = 'opto-off,opto-on,after-opto'; % columns: opto-off, opto-on, after-opto


% CONTROL PROTOCOL: OPTO-OFF
% search for max (within the PF) and min values (outside PF) and their positions in the position tuning curve
for j = 1:1:length(PlaceCells)
    i = PlaceCells(j); 
    PosTuningCirc = [sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig(i,1:nBins) sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig(i,1:nBins)];
    PFInAct = NaN(1,nBins);
    PFOutAct = NaN(1,nBins);
    for n = 1:1:size(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,:),2)
        if isnan(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,n)) 
           continue 
        end
        PFStartBin = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,n);
        PFEndBin = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldStartBin(i,n)+sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceFieldBinLength(i,n);
        PFInAct(PFStartBin:PFEndBin) = PosTuningCirc(PFStartBin:PFEndBin);
    end
    PFOutAct(isnan(PFInAct)) = PosTuningCirc(isnan(PFInAct));
    PFOutActScheme(i,1:nBins) = ~isnan(PFOutAct);
    % If the peak was in the discarded first part, discard ROIs data
    if isnan(max(PFInAct(1:nBins-1)))
        continue
    elseif find(PFInAct==(max(PFInAct(1:nBins-1)))) < DiscardBins + 1
        continue
    end
    % Discard the beginning of data because of optical stimulation delay    
    PFInAct(1:DiscardBins) = NaN;
    PFOutAct(1:DiscardBins) = NaN;
    % search for peak in PF and minimum (baseline) outside PF
    PeakPFin(i,1) = max(PFInAct(1:nBins-1)); % does not search in last bin, because it is problematic if highest peak is at last bin and I want to pool data with neighboring bins. In this case it will take the bin before or after
    PeakPFinPos(i,1) = find(PFInAct==PeakPFin(i,1),1);
    BaselinePFout(i,1) = min(PFOutAct(1:nBins-1));
    BaselinePFoutPos(i,1) = find(PFOutAct==BaselinePFout(i,1),1);
    % take the data of the two neighboring bins and average all three bins
    PeakPFin3BinsMean(i,1) = nanmean(PosTuningCirc(PeakPFinPos(i,1)-1:PeakPFinPos(i,1)+1));
    BaselinePFout3BinsMean(i,1) = nanmean(PosTuningCirc(BaselinePFoutPos(i,1)-1:BaselinePFoutPos(i,1)+1));    
    % Another Baseline calculation: take the 3 precentile value of the postition tuning curve 
    BaselinePFout3percentile(i,1) = prctile(PFOutAct,3);
        
end


% OPTO-ON PROTOCOL
% given the position of the peak in PF, and baseline, I calculate the activity values. 
%(IT is not done: I check if the peak of PF shifted a few bins. If yes, I will use the new position)

for j = 1:1:length(PlaceCells)
    i = PlaceCells(j); 
    % If the peak was in the discarded first part, discard ROIs data
    if isnan(PeakPFinPos(i,1))
        continue
    end
    PFOutAct = NaN(1,nBins);
    PosTuningCirc = [sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig(i,1:nBins) sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig(i,1:nBins)];
    PeakPFinPosCtr = PeakPFinPos(i,1);
    BaselinePFoutPosCtr = BaselinePFoutPos(i,1);
    
    PeakPFin(i,2) = PosTuningCirc(PeakPFinPosCtr); % does not search in last bin, because it is problematic if highest peak is at last bin and I want to pool data with neighboring bins. In this case it will take the bin before or after
    BaselinePFout(i,2) = PosTuningCirc(BaselinePFoutPosCtr);
    % take the data of the two neighboring bins and average all three bins
    PeakPFin3BinsMean(i,2) = nanmean(PosTuningCirc(PeakPFinPos(i,1)-1:PeakPFinPos(i,1)+1));
    BaselinePFout3BinsMean(i,2) = nanmean(PosTuningCirc(BaselinePFoutPos(i,1)-1:BaselinePFoutPos(i,1)+1)); 
    PosT = PosTuningCirc(1:nBins);
    PFOutAct(PFOutActScheme(i,:)==1) = PosT(PFOutActScheme(i,:)==1);
    BaselinePFout3percentile(i,2) = prctile(PFOutAct,3);

    %{
    % plot position tuning curve with max peak in PF and min in baseline
    figure('Color','white','visible','on')
    hold on
    % Opto-off
    plot(sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig(i,1:nBins),'-','Color','blue')
    plot(PeakPFinPos(i,1),PeakPFin(i,1),'o','Color','red');
    plot(PeakPFinPos(i,1),PeakPFin3BinsMean(i,1),'o','Color',[0.7, 0.4, 0.2]);
    plot(BaselinePFoutPos(i,1),BaselinePFout(i,1),'o','Color','blue');
    plot(BaselinePFoutPos(i,1),BaselinePFout3BinsMean(i,1),'o','Color','black');
    plot(0,BaselinePFout3BinsMean(i,1),'o','Color',[0, 0.5, 0]);
    % Opto-on
    plot(sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig(i,1:nBins),'-','Color','red')
    plot(PeakPFinPos(i,1),PeakPFin(i,2),'x','Color','red');
    plot(PeakPFinPos(i,1),PeakPFin3BinsMean(i,2),'x','Color',[0.7, 0.4, 0.2]);
    plot(BaselinePFoutPos(i,1),BaselinePFout(i,2),'x','Color','blue');
    plot(BaselinePFoutPos(i,1),BaselinePFout3BinsMean(i,2),'x','Color','black');
    plot(0,BaselinePFout3BinsMean(i,2),'x','Color',[0, 0.5, 0]);
    %axis([0 Max*1.1 0 1.1*Max]); 
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(i)));
    xlabel('Position on the wheel (cm)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('Activity (dF/F)');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(i),'PeakInPF_BaselineOutPF');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
%close all


% AFTER-OPTO PROTOCOL
% given the position of the peak in PF, and baseline, I calculate the activity values. 

for j = 1:1:length(PlaceCells)
    i = PlaceCells(j); 
    % If the peak was in the discarded first part, discard ROIs data
    if isnan(PeakPFinPos(i,1))
        continue
    end
    PosTuningCirc = [sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig(i,1:nBins) sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig(i,1:nBins)];
    PeakPFinPosCtr = PeakPFinPos(i,1);
    BaselinePFoutPosCtr = BaselinePFoutPos(i,1);
    
    PeakPFin(i,3) = PosTuningCirc(PeakPFinPosCtr); % does not search in last bin, because it is problematic if highest peak is at last bin and I want to pool data with neighboring bins. In this case it will take the bin before or after
    BaselinePFout(i,3) = PosTuningCirc(BaselinePFoutPosCtr);
    % take the data of the two neighboring bins and average all three bins
    PeakPFin3BinsMean(i,3) = nanmean(PosTuningCirc(PeakPFinPos(i,1)-1:PeakPFinPos(i,1)+1));
    BaselinePFout3BinsMean(i,3) = nanmean(PosTuningCirc(BaselinePFoutPos(i,1)-1:BaselinePFoutPos(i,1)+1));    
    PosT = PosTuningCirc(1:nBins);
    PFOutAct(PFOutActScheme(i,:)==1) = PosT(PFOutActScheme(i,:)==1);
    BaselinePFout3percentile(i,3) = prctile(PFOutAct,3);
    
end

sData.gainMod_RegressionForTwoPoints = struct;
sData.gainMod_RegressionForTwoPoints.DiscardedBins = DiscardCmBeginning;
sData.gainMod_RegressionForTwoPoints.PlaceCellsUsed = find(~isnan(PeakPFin(:,1)));
sData.gainMod_RegressionForTwoPoints.PeakPFin = PeakPFin;
sData.gainMod_RegressionForTwoPoints.PeakPFin3BinsMean = PeakPFin3BinsMean; % columns: opto-off, opto-on, after-opto
sData.gainMod_RegressionForTwoPoints.PeakPFinPos = PeakPFinPos;
sData.gainMod_RegressionForTwoPoints.BaselinePFout3BinsMean = BaselinePFout3BinsMean; 
sData.gainMod_RegressionForTwoPoints.BaselinePFoutPos = BaselinePFoutPos;
sData.gainMod_RegressionForTwoPoints.BaselinePFout3percentile = BaselinePFout3percentile;
sData.gainMod_RegressionForTwoPoints.LabelsForColumns = 'opto-off,opto-on,after-opto'; % columns: opto-off, opto-on, after-opto

sData.gainMod_RegressionForTwoPoints.means_OffOnAfter = struct;
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin(1,1) = nanmean(PeakPFin(:,1));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin(1,2) = nanmean(PeakPFin(:,2));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin(1,3) = nanmean(PeakPFin(:,3));

sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin3BinsMean(1,1) = nanmean(PeakPFin3BinsMean(:,1));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin3BinsMean(1,2) = nanmean(PeakPFin3BinsMean(:,2));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.PeakPFin3BinsMean(1,3) = nanmean(PeakPFin3BinsMean(:,3));

sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3BinsMean(1,1) = nanmean(BaselinePFout3BinsMean(:,1));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3BinsMean(1,2) = nanmean(BaselinePFout3BinsMean(:,2));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3BinsMean(1,3) = nanmean(BaselinePFout3BinsMean(:,3));

sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3percentile(1,1) = nanmean(BaselinePFout3percentile(:,1));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3percentile(1,2) = nanmean(BaselinePFout3percentile(:,2));
sData.gainMod_RegressionForTwoPoints.means_OffOnAfter.BaselinePFout3percentile(1,3) = nanmean(BaselinePFout3percentile(:,3));

% saving
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


end