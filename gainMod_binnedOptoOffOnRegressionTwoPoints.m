function sData = gainMod_binnedOptoOffOnRegressionTwoPoints(sData,savePath)

%savePath = 'E:\ANALYSIS\m8061-VIPArch-Thy1-GCamp\m8061-20200523-00-second-anal\Imaging\GainModTwoPointRegressionFits';

%mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'GainMod-RegressionTwoPoints');
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\GainMod-RegressionTwoPoints');
DiscardBins = sData.gainMod_RegressionForTwoPoints.DiscardedBins;
nROIs = sData.imdata.nROIs;

% linear fit for place cell rois on the previous plot
sData.sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.note = 'Equation for each place ROI: binned activity pairs (opto-off vs opto-on) are on the plot, one for peak in PF one for baseline, equation steepness, shift, regression and p value are determined.';
nPCRoi = size(sData.gainMod_RegressionForTwoPoints.PlaceCellsUsed,1);
PlaceCells = sData.gainMod_RegressionForTwoPoints.PlaceCellsUsed;
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.nPCRoi = nPCRoi;
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqSteepness = NaN(nROIs,1); % steepness of equation
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqShift = NaN(nROIs,1); % shift of equation
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.Regression2 = NaN(nROIs,1);
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.SignificanceP = NaN(nROIs,1);

% calculate linear fit for data and plot the opto-off vs opto-on amplitude value for each roi for each bin
for i = 1:1:nPCRoi
    roi = PlaceCells(i);
    PosTunOff = [sData.gainMod_RegressionForTwoPoints.BaselinePFout3BinsMean(roi,1) sData.gainMod_RegressionForTwoPoints.PeakPFin3BinsMean(roi,1)];
    PosTunOn = [sData.gainMod_RegressionForTwoPoints.BaselinePFout3BinsMean(roi,2) sData.gainMod_RegressionForTwoPoints.PeakPFin3BinsMean(roi,2)];
    Max = max(max(PosTunOff),max(PosTunOn));
    [fitX, fitY, R2, P, a, b] = linFit(PosTunOff,PosTunOn);
    sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.Regression2(roi) = R2;
    sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.SignificanceP(roi) = P;
    sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqSteepness(roi) = a;
    sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqShift(roi) = b;

    %{
    figure('Color','white','visible','on')
    hold on
    plot(PosTunOff,PosTunOn,'o')
    plot(0:1,0:1)
    plot(fitX,fitY,'--','Color','black');
    axis([0 Max*1.1 0 1.1*Max]); 
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Mean binned activity in Opto-Off trials (\DeltaF/F)'); ax = gca; ax.TickDir = 'out'; 
    ylabel('Mean binned activity in Opto-On trials (\DeltaF/F)');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-BinnedOffVsOn');
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

    if rem(i,20)==0
       close all 
    end
    %}
end

sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqSteepness = nanmean(sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqSteepness);
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqShift = nanmean(sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.EqShift);
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqRegression2 = nanmean(sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.Regression2);
sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqSignificanceP = nanmean(sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.SignificanceP);

%plot a mean regression figure for the session
x = (0:0.1:1);
m = sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqSteepness;
b = sData.gainMod_RegressionForTwoPoints.binnedOptoOffOnRegression.meanEqShift;
y = m*x+b;
% plot normal figure
figure('Color','white','Position',[100 100 350 300])
hold on
plot(0:1,0:1,'Color','red','LineWidth',0.5)
plot(x,y,'--','Color','black','LineWidth',1); % plot mean regression line
axis([0 1 0 1]); 
xlabel('Opto-Off'); 
ylabel('Opto-On');
ax = gca; ax.TickDir = 'out'; ax.FontSize = 12; 
xticks([0 0.2 0.4 0.6 0.8 1]); yticks([0 0.2 0.4 0.6 0.8 1]);
fname = strcat(sData.sessionInfo.fileID,'-MeanRegression');
exportgraphics(gcf,fullfile(savePath,[fname '.png']),'Resolution',300);
exportgraphics(gcf,fullfile(savePath,[fname '.pdf']),'Resolution',300);
saveas(gcf,(fullfile(savePath,[fname '.pdf'])));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% plot mini figure
figure('Color','white','Position',[100 100 200 200])
hold on
plot(0:1,0:1,'Color','red','LineWidth',0.5)
plot(x,y,'--','Color','black','LineWidth',1); % plot mean regression line
axis([0 1 0 1]); 
xlabel('Opto-Off'); 
ylabel('Opto-On');
ax = gca; ax.TickDir = 'out'; ax.FontSize = 16; 
xticks([0 1]); yticks([0 1]);
fname = strcat(sData.sessionInfo.fileID,'-MeanRegressionMini');
exportgraphics(gcf,fullfile(savePath,[fname '.png']),'Resolution',300);
exportgraphics(gcf,fullfile(savePath,[fname '.pdf']),'Resolution',300);
saveas(gcf,(fullfile(savePath,[fname '.pdf'])));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% Save file to same path where other files can be found 
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end