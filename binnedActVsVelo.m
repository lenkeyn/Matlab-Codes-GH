function sData = binnedActVsVelo(sData,binToNorm) % bin to norm means that normalize Y value to that bins value for both velocity and ROI act

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'VeloNormAct');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\VeloNormAct'); 
% binToNorm = 20;

% cue positions: %{
C1A = 23; C1B = 29; % velcro
C2A = 43; C2B = 49; % hot glue
C3A = 63; C3B = 69; % hot glue
C4A = 83; C4B = 89; % velcro

% COMPARE binned ROI activity with mean velocity curve
% normalized velocity
MaxVelo = max(sData.behavior.binning.meanVeloBinned,[],2);
MinVelo = min(sData.behavior.binning.meanVeloBinned,[],2);
NormVeloBinned = (sData.behavior.binning.meanVeloBinned - MinVelo) ./  (MaxVelo - MinVelo); 

% normalization of MeanBinnedActivity plots for  ROIs 
MaxAct = max(sData.imdata.binned.MeanRoiAct,[],2);
MinAct = min(sData.imdata.binned.MeanRoiAct,[],2);
sData.imdata.binned.analysis.MeanRoiActNorm = (sData.imdata.binned.MeanRoiAct - MinAct) ./  (MaxAct - MinAct); 
MeanMeanRoiAct = mean(sData.imdata.binned.analysis.MeanRoiActNorm,1);

% For Grand Average I normalize activity max at 20th bin
%MaxAct = max(MeanMeanRoiAct,[],2);
MinAct = min(MeanMeanRoiAct,[],2);
MaxAct = MeanMeanRoiAct(binToNorm);
sData.imdata.binned.analysis.MeanMeanRoiActNormAll = (MeanMeanRoiAct - MinAct) ./  (MaxAct - MinAct); 

% compare running profile and average imaging of all ROIs
figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
plot(Xaxis,NormVeloBinned(1:sData.behavior.meta.nBins)); hold on
plot(Xaxis,sData.imdata.binned.analysis.MeanMeanRoiActNormAll(1:sData.behavior.meta.nBins)); hold on
ylim = get(gca,'YLim')-0.01;
line([C1A C1A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C1B C1B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
line([C2A C2A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C2B C2B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C3A C3A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C3B C3B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C4A C4A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C4B C4B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
title(strcat(sData.sessionInfo.fileID,' Velocity vs Mean Activity (all ROIs)'));
legend('Norm. velocity','Norm. activity of all ROIs','Location','south');
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Velocity / Mean activity');
fname = strcat(sData.sessionInfo.fileID,'-VelocityVsMeanActivity-allROIs');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% collect neuronal soma ROIs into array
somaArray = NaN(sData.imdata.nROIs,1);
for i = 1:1:sData.imdata.nROIs
    if strcmp(sData.imdata.roiArray(i).group,'Neuronal Soma')
        somaArray(i) = i;
    end
end
MeanMeanRoiActNormSoma = nanmean(sData.imdata.binned.analysis.MeanRoiActNorm(somaArray(~isnan(somaArray)),:));
MaxAct = MeanMeanRoiActNormSoma(binToNorm);
MinAct = min(MeanMeanRoiActNormSoma,[],2);
sData.imdata.binned.analysis.MeanMeanRoiActNormSoma = (MeanMeanRoiActNormSoma - MinAct) ./  (MaxAct - MinAct);
% compare running profile and average imaging
figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
plot(Xaxis,NormVeloBinned(1:sData.behavior.meta.nBins)); hold on
plot(Xaxis,sData.imdata.binned.analysis.MeanMeanRoiActNormSoma(1:sData.behavior.meta.nBins)); hold on
ylim = get(gca,'YLim')-0.01;
line([C1A C1A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C1B C1B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
line([C2A C2A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C2B C2B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C3A C3A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C3B C3B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C4A C4A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C4B C4B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
title(strcat(sData.sessionInfo.fileID,' Velocity vs Mean Activity (soma)'));
legend('Norm. velocity','Norm. activity of somatic ROIs','Location','south');
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Velocity / Mean activity');
fname = strcat(sData.sessionInfo.fileID,'-VelocityVsMeanActivity-somaROIs');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

% collect neuronal dendritic ROIs into array
% binToNorm = 10;
dendriteArray = NaN(sData.imdata.nROIs,1);
for i = 1:1:sData.imdata.nROIs
    if strcmp(sData.imdata.roiArray(i).group,'Neuronal Dendrite')
        dendriteArray(i) = i;
    end
end
MeanMeanRoiActNormDendrite = nanmean(sData.imdata.binned.analysis.MeanRoiActNorm(dendriteArray(~isnan(dendriteArray)),:));
MaxAct = MeanMeanRoiActNormDendrite(binToNorm);
MinAct = min(MeanMeanRoiActNormDendrite,[],2);
sData.imdata.binned.analysis.MeanMeanRoiActNormDendrite = (MeanMeanRoiActNormDendrite - MinAct) ./  (MaxAct - MinAct);
% compare running profile and average imaging
figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
plot(Xaxis,NormVeloBinned(1:sData.behavior.meta.nBins)); hold on
plot(Xaxis,sData.imdata.binned.analysis.MeanMeanRoiActNormDendrite(1:sData.behavior.meta.nBins)); hold on
ylim = get(gca,'YLim')-0.01;
line([C1A C1A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C1B C1B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
line([C2A C2A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C2B C2B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C3A C3A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C3B C3B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C4A C4A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C4B C4B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
title(strcat(sData.sessionInfo.fileID,' Velocity vs Mean Activity (dendrites)'));
legend('Norm. velocity','Norm. activity of dendrite ROIs','Location','south');
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Velocity / Mean activity');
fname = strcat(sData.sessionInfo.fileID,'-VelocityVsMeanActivity-dendrites');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% collect neuropil ROIs into array
neuropilArray = NaN(sData.imdata.nROIs,1);
for i = 1:1:sData.imdata.nROIs
    if strcmp(sData.imdata.roiArray(i).group,'Neuropil')
        neuropilArray(i) = i;
    end
end
MeanMeanRoiActNormNeuropil = nanmean(sData.imdata.binned.analysis.MeanRoiActNorm(neuropilArray(~isnan(neuropilArray)),:));
MaxAct = MeanMeanRoiActNormNeuropil(binToNorm);
MinAct = min(MeanMeanRoiActNormNeuropil,[],2);
sData.imdata.binned.analysis.MeanMeanRoiActNormNeuropil = (MeanMeanRoiActNormNeuropil - MinAct) ./  (MaxAct - MinAct);
% compare running profile and average imaging
figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
plot(Xaxis,NormVeloBinned(1:sData.behavior.meta.nBins)); hold on
plot(Xaxis,sData.imdata.binned.analysis.MeanMeanRoiActNormNeuropil(1:sData.behavior.meta.nBins)); hold on
ylim = get(gca,'YLim')-0.01;
line([C1A C1A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C1B C1B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
line([C2A C2A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C2B C2B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C3A C3A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C3B C3B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C4A C4A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C4B C4B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
title(strcat(sData.sessionInfo.fileID,' Velocity vs Mean Activity (neuropil)'));
legend('Norm. velocity','Norm. activity of neuropil ROIs','Location','south');
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Velocity / Mean activity');
fname = strcat(sData.sessionInfo.fileID,'-VelocityVsMeanActivity-Neuropil');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% collect gliopil ROIs into array (strange spiny dendrite, seems PC dendrite)
% binToNorm = 60;
strangeDendriteArray = NaN(sData.imdata.nROIs,1);
for i = 1:1:sData.imdata.nROIs
    if strcmp(sData.imdata.roiArray(i).group,'Gliopill')
        strangeDendriteArray(i) = i;
    end
end
MeanMeanRoiActNormstrangeDendrite = nanmean(sData.imdata.binned.analysis.MeanRoiActNorm(strangeDendriteArray(~isnan(strangeDendriteArray)),:));
MaxAct = MeanMeanRoiActNormstrangeDendrite(binToNorm);
MinAct = min(MeanMeanRoiActNormstrangeDendrite,[],2);
sData.imdata.binned.analysis.MeanMeanRoiActNormstrangeDendrite = (MeanMeanRoiActNormstrangeDendrite - MinAct) ./  (MaxAct - MinAct);
% compare running profile and average imaging
figure();
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
plot(Xaxis,NormVeloBinned(1:sData.behavior.meta.nBins)); hold on
plot(Xaxis,sData.imdata.binned.analysis.MeanMeanRoiActNormstrangeDendrite(1:sData.behavior.meta.nBins)); hold on
ylim = get(gca,'YLim')-0.01;
line([C1A C1A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C1B C1B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
line([C2A C2A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C2B C2B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C3A C3A],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on; line([C3B C3B],[0 ylim(2)],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',0.5); hold on;
line([C4A C4A],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on; line([C4B C4B],[0 ylim(2)],'Color','black','LineStyle','--','LineWidth',0.5); hold on;
title(strcat(sData.sessionInfo.fileID,' Velocity vs Mean Activity (spiny dendrite)'));
legend('Norm. velocity','Norm. activity of spiny dendrite ROIs','Location','north');
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Velocity / Mean activity');
fname = strcat(sData.sessionInfo.fileID,'-VelocityVsMeanActivity-SpinyDendrite');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end
