function sData = gainMod_binnedOptoOffOnRegressionPCinAnyProtocols(sData,figGeneration,DiscardCmBeginning)

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'GainModPCinAnyProt');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\GainModPCinAnyProt');
DiscardBins = round(DiscardCmBeginning/sData.behavior.meta.binSize);
nROIs = sData.imdata.nROIs;
%nBins = sData.behavior.meta.nBins - DiscardBins; % discard the first 20 cm, since residual dFF can generate false transients
SteepnessCtrSession = 0.658; % compare opto sessions steepness to ctr

% linear fit for place cell rois on the previous plot
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.note = 'Equation for each place ROI: binned activity pairs (opto-off vs opto-on) are on the plot, equation steepness, shift, regression and p value are determined.';
nPCRoi = size(sData.imdata.MaoPC_Opto_dff.PlaceCellinAnyProt,1);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.nPCRoi = nPCRoi;
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqSteepness = NaN(nROIs,1); % steepness of equation
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqShift = NaN(nROIs,1); % shift of equation
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.Regression2 = NaN(nROIs,1);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.SignificanceP = NaN(nROIs,1);

% calculate linear fit for data and plot the opto-off vs opto-on amplitude value for each roi for each bin
for i = 1:1:nPCRoi
    roi = sData.imdata.MaoPC_Opto_dff.PlaceCellinAnyProt(i);
    PosTunOff = sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuning(roi,DiscardBins+1:sData.behavior.meta.nBins);
    PosTunOn = sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuning(roi,DiscardBins+1:sData.behavior.meta.nBins);
    Max = max(max(PosTunOff),max(PosTunOn));
    [fitX, fitY, R2, P, a, b] = linFit(PosTunOff,PosTunOn);
    sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.Regression2(roi) = R2;
    sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.SignificanceP(roi) = P;
    sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqSteepness(roi) = a;
    sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqShift(roi) = b;

    if figGeneration > 0
        figure('Color','white','visible','off')
        hold on
        plot(PosTunOff,PosTunOn,'o')
        plot(0:1,0:1)
        plot(fitX,fitY,'--','Color','black');
        axis([0 Max*1.1 0 1.1*Max]); 
        title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
        xlabel('Mean binned activity in Opto-Off trials (dF/F)'); ax = gca; ax.TickDir = 'out'; 
        ylabel('Mean binned activity in Opto-On trials (dF/F)');
        fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'BinnedOffVsOn');
        savefig(fullfile(savePath,fname));
        saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

        if rem(i,20)==0
           close all 
        end
    end
    
end

sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.meanEqSteepness = nanmean(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqSteepness);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.meanEqShift = nanmean(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqShift);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.meanEqRegression2 = nanmean(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.Regression2);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.meanEqSignificanceP = nanmean(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.SignificanceP);
sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.ratioEqSteepnessAboveCtr = length(find(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqSteepness>SteepnessCtrSession))/length(find(sData.gainModulationPCinAnyProtocol.binnedOptoOffOnRegression.EqSteepness<SteepnessCtrSession));


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end