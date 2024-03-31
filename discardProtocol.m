function sData = discardProtocol(sData,discardProt)

% selecting only one protocol if there are more:
 % keepProt = 2; discardProt = 3; % #1 is control,#2-3-4... is the optically stimulated

sData.behavior.opto.OptoOnTrials(sData.behavior.optoMoreProts.OptoStimProtTrialsReal == discardProt) = NaN;
sData.behavior.opto.OptoOnTrialsIndices = find(sData.behavior.opto.OptoOnTrials == 1);
sData.behavior.opto.AfterOptoTrials(sData.behavior.opto.OptoOnTrialsIndices+1) = NaN;
sData.behavior.opto.AfterOptoTrialsIndices = find(sData.behavior.opto.AfterOptoTrials == 1);


%%% Sort the licking data
sData.behavior.opto.LickMatrix_OptoOff = sData.behavior.binning.lickPerCmBinned(sData.behavior.opto.OptoOffTrialsIndices,:);
sData.behavior.opto.LickMatrix_OptoOn = sData.behavior.binning.lickPerCmBinned(sData.behavior.opto.OptoOnTrialsIndices,:);
sData.behavior.opto.LickMatrix_AfterOpto = sData.behavior.binning.lickPerCmBinned(sData.behavior.opto.AfterOptoTrialsIndices,:);

sData.behavior.opto.MeanLick_OptoOff = nanmean(sData.behavior.opto.LickMatrix_OptoOff,1);
sData.behavior.opto.MeanLick_OptoOn = nanmean(sData.behavior.opto.LickMatrix_OptoOn,1);
sData.behavior.opto.MeanLick_OptoAfter = nanmean(sData.behavior.opto.LickMatrix_AfterOpto,1);

%%%% Sort speed data
sData.behavior.opto.VeloMatrix_OptoOff = sData.behavior.binning.veloBinned(sData.behavior.opto.OptoOffTrialsIndices,:);
sData.behavior.opto.VeloMatrix_OptoOn = sData.behavior.binning.veloBinned(sData.behavior.opto.OptoOnTrialsIndices,:);
sData.behavior.opto.VeloMatrix_AfterOpto = sData.behavior.binning.veloBinned(sData.behavior.opto.AfterOptoTrialsIndices,:);

sData.behavior.opto.MeanVelo_OptoOff = mean(sData.behavior.opto.VeloMatrix_OptoOff,1);
sData.behavior.opto.MeanVelo_OptoOn = mean(sData.behavior.opto.VeloMatrix_OptoOn,1);
sData.behavior.opto.MeanVelo_OptoAfter = mean(sData.behavior.opto.VeloMatrix_AfterOpto,1);

%%% HitRate
PosNaN = find(isnan(sData.behavior.details.hitTrials));
sData.behavior.opto.ActiveHitRate_OptoOff = sum(sData.behavior.details.hitTrials(sData.behavior.opto.OptoOffTrialsIndices(sData.behavior.opto.OptoOffTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(sData.behavior.opto.OptoOffTrialsIndices,1))));
sData.behavior.opto.ActiveHitRate_OptoOn = sum(sData.behavior.details.hitTrials(sData.behavior.opto.OptoOnTrialsIndices(sData.behavior.opto.OptoOnTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(sData.behavior.opto.OptoOnTrialsIndices,1))));
sData.behavior.opto.ActiveHitRate_OptoAfter = sum(sData.behavior.details.hitTrials(sData.behavior.opto.AfterOptoTrialsIndices(sData.behavior.opto.AfterOptoTrialsIndices<PosNaN),1)/nnz(~isnan(sData.behavior.details.hitTrials(sData.behavior.opto.AfterOptoTrialsIndices,1))));

end