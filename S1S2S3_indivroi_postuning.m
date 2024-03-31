% mean pos tuning in one figure if S1, S2, S3. Open sData for all session, and write Roi ID
    roi = 100;
    savePath = 'C:\MATLAB\SAVE';
    nBins = sData1.behavior.meta.nBins;
    ROI(1,:)= sData1.imdata.binned.MeanRoiAct(roi,1:nBins);
    ROI(2,:)= sData2.imdata.binned.MeanRoiAct(roi,1:nBins);
    ROI(3,:)= sData3.imdata.binned.MeanRoiAct(roi,1:nBins);
    
    figure();
    Xaxis = sData1.behavior.meta.binSize:sData1.behavior.meta.binSize:sData1.behavior.meta.binSize*nBins;
    Ymax = 0.3; %(max(nanmean(MeanRoiActDenoised)))*1.1;
    Ymin = -0.05; %(min(nanmean(MeanRoiActDenoised)))*0.9;
    plot(Xaxis,ROI(1,1:nBins),'LineWidth',2); hold on;
    plot(Xaxis,ROI(2,1:nBins),'LineWidth',2);
    plot(Xaxis,ROI(3,1:nBins),'LineWidth',2);
    maxpos1 = find(ROI(1,1:nBins)== max(ROI(1,1:nBins)))*sData1.behavior.meta.binSize;
    maxpos2 = find(ROI(2,1:nBins)== max(ROI(2,1:nBins)))*sData1.behavior.meta.binSize;
    maxpos3 = find(ROI(3,1:nBins)== max(ROI(3,1:nBins)))*sData1.behavior.meta.binSize;
    line([maxpos1 maxpos1],[Ymin max(ROI(1,1:nBins))],'Color',[0, 0.4470, 0.7410],'LineStyle','--','LineWidth',2); 
    line([maxpos2 maxpos2],[Ymin max(ROI(2,1:nBins))],'Color',[0.8500, 0.3250, 0.0980],'LineStyle','--','LineWidth',2); 
    line([maxpos3 maxpos3],[Ymin max(ROI(3,1:nBins))],'Color',[0.9290, 0.6940, 0.1250],'LineStyle','--','LineWidth',2); 
    %line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 Ymin Ymax]); % ceil(Ymax)
    title(strcat(sData1.sessionInfo.fileID,'-ROI-',num2str(roi),'-Position tuning'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    legend('RF-1','GOL','RF-2','location','north');
    fname = strcat(sData1.sessionInfo.fileID,'ROI',num2str(roi),'Pos-tuning-S1-S2-S3');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));