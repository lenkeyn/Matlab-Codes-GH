% ROUGH ANALYSIS ON THE CELLULAR EFFET
%CALCUALTE THE MEAN OF THE POSITION TUNING CURVE AND COMPARE THIS VALUE IN OPTO-OFF AND ON CONDITIONS FOR ALL CELLS AND PLACE CELLS

% ONE OPTO - ALL CELLS  
MeanPosTuningOff = nanmean(sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig,1);
MeanPosTuningOpto = nanmean(sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig,1);
MeanPosTuningAfter = nanmean(sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig,1);

MMAllCellsPosTuningFull(1,1) = mean(MeanPosTuningOff);
MMAllCellsPosTuningFull(1,2) = mean(MeanPosTuningOpto);
MMAllCellsPosTuningFull(1,3) = mean(MeanPosTuningAfter);

MMAllCellsPosTuningBin680(1,1) = mean(MeanPosTuningOff(:,11:end));
MMAllCellsPosTuningBin680(1,2) = mean(MeanPosTuningOpto(:,11:end));
MMAllCellsPosTuningBin680(1,3) = mean(MeanPosTuningAfter(:,11:end));

% ONE OPTO - ONLY PLACE CELLS
MeanPosTuningPCOff = mean(sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuningOrig(sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells,:),1);
MeanPosTuningPCOpto = mean(sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuningOrig(sData.imdata.MaoPC_Opto_dff.OptoOn.PlaceCells,:),1);
MeanPosTuningPCAfter = mean(sData.imdata.MaoPC_Opto_dff.OptoAfter.PosTuningOrig(sData.imdata.MaoPC_Opto_dff.OptoAfter.PlaceCells,:),1);

MMPCPosTuningFull(1,1) = mean(MeanPosTuningPCOff);
MMPCPosTuningFull(1,2) = mean(MeanPosTuningPCOpto);
MMPCPosTuningFull(1,3) = mean(MeanPosTuningPCAfter);

MMPCPosTuningBin680(1,1) = mean(MeanPosTuningOff(:,11:end));
MMPCPosTuningBin680(1,2) = mean(MeanPosTuningOpto(:,11:end));
MMPCPosTuningBin680(1,3) = mean(MeanPosTuningAfter(:,11:end));


% MORE OPTO ALL CELLS
MeanPosTuningOff = nanmean(sData.imdata.MaoOptoMoreProt_dff{1, 1}.PosTuningOrig,1);
MeanPosTuningOpto1 = nanmean(sData.imdata.MaoOptoMoreProt_dff{1, 2}.PosTuningOrig,1);
MeanPosTuningOpto2 = nanmean(sData.imdata.MaoOptoMoreProt_dff{1, 3}.PosTuningOrig,1);
MeanPosTuningAfter = nanmean(sData.imdata.MaoOptoMoreProt_dff{1, 5}.PosTuningOrig,1);

MMAllCellsPosTuningFull(1,1) = mean(MeanPosTuningOff);
MMAllCellsPosTuningFull(1,2) = mean(MeanPosTuningOpto2);
MMAllCellsPosTuningFull(1,3) = mean(MeanPosTuningOpto1);
MMAllCellsPosTuningFull(1,4) = mean(MeanPosTuningAfter);

MMAllCellsPosTuningBin680(1,1) = mean(MeanPosTuningOff(:,11:end));
MMAllCellsPosTuningBin680(1,2) = mean(MeanPosTuningOpto2(:,11:end));
MMAllCellsPosTuningBin680(1,3) = mean(MeanPosTuningOpto1(:,11:end));
MMAllCellsPosTuningBin680(1,4) = mean(MeanPosTuningAfter(:,11:end));

% MORE OPTO ONLY PLACE CELLS
MeanPosTuningPCOff = mean(sData.imdata.MaoOptoMoreProt_dff{1, 1}.PosTuningOrig(sData.imdata.MaoOptoMoreProt_dff{1, 1}.PlaceCells,:),1);
MeanPosTuningPCOpto1 = mean(sData.imdata.MaoOptoMoreProt_dff{1, 2}.PosTuningOrig(sData.imdata.MaoOptoMoreProt_dff{1, 2}.PlaceCells,:),1);
MeanPosTuningPCOpto2 = mean(sData.imdata.MaoOptoMoreProt_dff{1, 3}.PosTuningOrig(sData.imdata.MaoOptoMoreProt_dff{1, 3}.PlaceCells,:),1);
MeanPosTuningPCAfter = mean(sData.imdata.MaoOptoMoreProt_dff{1, 5}.PosTuningOrig(sData.imdata.MaoOptoMoreProt_dff{1, 5}.PlaceCells,:),1);

MMPCPosTuningFull(1,1) = mean(MeanPosTuningPCOff);
MMPCPosTuningFull(1,2) = mean(MeanPosTuningPCOpto2);
MMPCPosTuningFull(1,3) = mean(MeanPosTuningPCOpto1);
MMPCPosTuningFull(1,4) = mean(MeanPosTuningPCAfter);

MMPCPosTuningBin680(1,1) = mean(MeanPosTuningOff(:,11:end));
MMPCPosTuningBin680(1,2) = mean(MeanPosTuningOpto2(:,11:end));
MMPCPosTuningBin680(1,3) = mean(MeanPosTuningOpto1(:,11:end));
MMPCPosTuningBin680(1,4) = mean(MeanPosTuningAfter(:,11:end));

