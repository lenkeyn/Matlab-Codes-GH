
%sumTable = NaN(20,40);

summary = struct;

% recording ID
n = 1;
sessionID = 70712020052700; % m8058-20200527-01  => 80582020052701

summary.meanROIActBinPC.OptoOff = mean(sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuning(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,:));
summary.meanROIActBinPC.OptoOn = mean(sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuning(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,:));
summary.meanROIActBinPC.OptoAfter = mean(sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuning(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,:));

summary.meanROIActBinPC.AUC_OptoOff = trapz(summary.meanROIActBinPC.OptoOff);
summary.meanROIActBinPC.AUC_OptoOn = trapz(summary.meanROIActBinPC.OptoOn);
summary.meanROIActBinPC.AUC_AfterOpto = trapz(summary.meanROIActBinPC.OptoAfter);

summary.meanROIActBinAllCells.OptoOff = nanmean(sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuning);
summary.meanROIActBinAllCells.OptoOn = nanmean(sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuning);
summary.meanROIActBinAllCells.OptoAfter = nanmean(sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuning);

summary.meanROIActBinAllCells.AUC_OptoOff = trapz(summary.meanROIActBinAllCells.OptoOff);
summary.meanROIActBinAllCells.AUC_OptoOn = trapz(summary.meanROIActBinAllCells.OptoOn);
summary.meanROIActBinAllCells.AUC_AfterOpto = trapz(summary.meanROIActBinAllCells.OptoAfter);

sData.analysisSummary = summary;
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

%%% Write data into summary table
sumTable(n,1) = sessionID;
sumTable(n,2) = str2double(erase(sData.mouseInfo.dateOfBirth,"."));
sumTable(n,3) = sData.behavior.opto.ActiveHitRate_OptoOff;
sumTable(n,4) = sData.behavior.opto.ActiveHitRate_OptoOn;
sumTable(n,5) = sData.behavior.opto.ActiveHitRate_OptoAfter;
sumTable(n,6) = sData.behavior.opto.optoStimStart;
sumTable(n,7) = sData.behavior.opto.optoStimEnd;
sumTable(n,8) = sData.behavior.performance2lick.meanQ50; % licks were deleted from the first 20% of trials
sumTable(n,9) = mean(sData.behavior.performance.velo.veloMax);
sumTable(n,10) = sData.behavior.performance2lick.MeanFirstLickCm;
sumTable(n,11) = sum(sData.behavior.opto.FailedOptoTrials(~isnan(sData.behavior.opto.FailedOptoTrials))); % number of failed trials
sumTable(n,12) = prctile(sData.imdata.binned.MeanMotCorrBinned,90);
sumTable(n,13) = sData.imdata.nROIs;
sumTable(n,14) = sData.imdata.nSamples;
sumTable(n,15) = sData.behavior.wheelLapImaging;
sumTable(n,16) = sData.imdata.signalExtractionOptions.spk_SNR;
sumTable(n,17) = sData.imdata.signalExtractionOptions.lam_pr;
sumTable(n,18) = sData.imdata.roiStat.meanPeakDff;
sumTable(n,19) = sData.imdata.roiStat.meanSignalToNoise;
sumTable(n,20) = sData.imdata.roiStat.meanActivityLevel;
sumTable(n,21) = size(sData.imdata.MaoPC_Opto_dff.LandmarkCells.LandmarkCellinAnyProt,1);   
sumTable(n,22) = size(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.PlaceCellsWithOnePFinAnyProt,1);
sumTable(n,23) = sData.imdata.MaoPC_Opto_dff.OptoOff.placeROINu;
sumTable(n,24) = sData.imdata.MaoPC_Opto_dff.OptoOn.placeROINu;
sumTable(n,25) = sData.imdata.MaoPC_Opto_dff.OptoAfter.placeROINu;
sumTable(n,26) = length(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOff.LandmarkCellsList);
sumTable(n,27) = length(sData.imdata.MaoPC_Opto_dff.LandmarkCells.OptoOn.LandmarkCellsList);
sumTable(n,28) = length(sData.imdata.MaoPC_Opto_dff.LandmarkCells.AfterOpto.LandmarkCellsList);
sumTable(n,29) = length(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff.PlaceCellsWithOnePFList);
sumTable(n,30) = length(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOn.PlaceCellsWithOnePFList);
sumTable(n,31) = length(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoAfter.PlaceCellsWithOnePFList);

sumTable(n,32) = summary.meanROIActBinPC.AUC_OptoOff;
sumTable(n,33) = summary.meanROIActBinPC.AUC_OptoOn;
sumTable(n,34) = summary.meanROIActBinPC.AUC_AfterOpto;

sumTable(n,35) = summary.meanROIActBinAllCells.AUC_OptoOff;
sumTable(n,36) = summary.meanROIActBinAllCells.AUC_OptoOn;
sumTable(n,37) = summary.meanROIActBinAllCells.AUC_AfterOpto;

sumTable(n,38) = sData.behavior.details.MaxBackTRack;