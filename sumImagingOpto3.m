function sumTable = sumImagingOpto3(sData)

sumTable = NaN(1,18);


%%% Write data into summary table
% PF-InOutRatio
sumTable(1,1) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,1);
sumTable(1,2) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,2);
sumTable(1,3) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.meanInMeanOutRatio_OffOnAfter(1,3);
sumTable(1,4) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,1);
sumTable(1,5) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,2);
sumTable(1,6) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakInMeanOutRatio_OffOnAfter(1,3);
sumTable(1,7) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,1);
sumTable(1,8) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,2); 
sumTable(1,9) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.peakActInPlaceField_OffOnAfter(1,3);
sumTable(1,10) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,1);
sumTable(1,11) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,2); 
sumTable(1,12) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.BinPosPeakActInPlaceField_OffOnAfter(1,3);

% regression
sumTable(1,13) = sData.gainModulationPCinLandmarkAndOnePFPC.binnedOptoOffOnRegression.meanEqSteepness;
sumTable(1,14) = sData.gainModulationPCinLandmarkAndOnePFPC.binnedOptoOffOnRegression.meanEqShift;
sumTable(1,15) = sData.gainModulationPCinLandmarkAndOnePFPC.binnedOptoOffOnRegression.meanEqRegression2;
sumTable(1,16) = sData.gainModulationPCinLandmarkAndOnePFPC.binnedOptoOffOnRegression.meanEqSignificanceP;

sumTable(1,17) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.InOnPerInOff;
sumTable(1,18) = sData.gainModulationPCinLandmarkAndOnePFPC.InOutMeanActRatio.means.OutOnPerOutOff;

end