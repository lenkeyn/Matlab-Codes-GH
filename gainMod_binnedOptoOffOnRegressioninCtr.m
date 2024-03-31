function sData = gainMod_binnedOptoOffOnRegressioninCtr(sData,figGeneration,DiscardCmBeginning)

%mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'GainModPCinCtr');
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\GainModPCinCtr');
DiscardBins = round(DiscardCmBeginning/sData.behavior.meta.binSize);
nROIs = sData.imdata.nROIs;
PlaceCells = sData.imdata.MaoPC_Opto_dff.OptoOff.PlaceCells; 
%nBins = sData.behavior.meta.nBins - DiscardBins; % discard the first 20 cm, since residual dFF can generate false transients
%SteepnessCtrSession = 0.658; % compare opto sessions steepness to ctr

% linear fit for place cell rois on the previous plot
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.note = 'Equation for each place ROI: binned activity pairs (opto-off vs opto-on) are on the plot, equation steepness, shift, regression and p value are determined.';
nPCRoi = size(PlaceCells,1);
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.nPCRoi = nPCRoi;
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqSteepness = NaN(nROIs,1); % steepness of equation
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqShift = NaN(nROIs,1); % shift of equation
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.Regression2 = NaN(nROIs,1);
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.SignificanceP = NaN(nROIs,1);

% calculate linear fit for data and plot the opto-off vs opto-on amplitude value for each roi for each bin
for i = 1:1:nPCRoi
    roi = PlaceCells(i);
    PosTunOff = sData.imdata.MaoPC_Opto_dff.OptoOff.PosTuning(roi,DiscardBins+1:sData.behavior.meta.nBins);
    PosTunOn = sData.imdata.MaoPC_Opto_dff.OptoOn.PosTuning(roi,DiscardBins+1:sData.behavior.meta.nBins);
    Max = max(max(PosTunOff),max(PosTunOn));
    [fitX, fitY, R2, P, a, b] = linFit(PosTunOff,PosTunOn);
    sData.gainModulationPCinCtr.binnedOptoOffOnRegression.Regression2(roi) = R2;
    sData.gainModulationPCinCtr.binnedOptoOffOnRegression.SignificanceP(roi) = P;
    sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqSteepness(roi) = a;
    sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqShift(roi) = b;

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

sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqSteepness = nanmean(sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqSteepness);
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqShift = nanmean(sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqShift);
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqRegression2 = nanmean(sData.gainModulationPCinCtr.binnedOptoOffOnRegression.Regression2);
sData.gainModulationPCinCtr.binnedOptoOffOnRegression.meanEqSignificanceP = nanmean(sData.gainModulationPCinCtr.binnedOptoOffOnRegression.SignificanceP);
%sData.gainModulationPCinCtr.binnedOptoOffOnRegression.ratioEqSteepnessAboveCtr = length(find(sData.gainModulationPCinCtr.binnedOptoOffOnRegression.EqSteepness>SteepnessCtrSession))/length(find(sData.gainModulationPCininCtr.binnedOptoOffOnRegression.EqSteepness<SteepnessCtrSession));


% Save file to same path where other files can be found 
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end