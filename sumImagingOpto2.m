function sumTable2 = sumImagingOpto2(sData)

sumTable2 = NaN(1,18);


%%% Write data into summary table
% PF-InOutRatio
sumTable2(1,1) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1);
sumTable2(1,2) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2);
sumTable2(1,3) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3);
sumTable2(1,4) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1);
sumTable2(1,5) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2);
sumTable2(1,6) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3);
sumTable2(1,7) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1);
sumTable2(1,8) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2); 
sumTable2(1,9) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3);
sumTable2(1,10) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1);
sumTable2(1,11) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2); 
sumTable2(1,12) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3);

% regression
sumTable2(1,13) = sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqSteepness;
sumTable2(1,14) = sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqShift;
sumTable2(1,15) = sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqRegression2;
sumTable2(1,16) = sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqSignificanceP;

sumTable2(1,17) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.InOnPerInOff;
sumTable2(1,18) = sData.gainModulationPCinCtr.InOutMeanActRatio.means.OutOnPerOutOff;


end