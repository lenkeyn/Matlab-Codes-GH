function sData  = CaDataOptoSortingdFFsData2(sData,FigVisible) 

% type: which data type to use, it can be 0 = RoidFF or 1 = RoidFF-SR (slow transients removed), 2= RoidFF-SR-LP (slow transients removed and transients are low pass filtered)
% IsOpto  = 1 : it was an optically stimulated session, IsOpto = 0 : sessions without stimulation

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct-Opto');
SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-Opto');
BinSize = sData.behavior.meta.binSize;
BinNu = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;

%Sorting binned trials for plotting
sData.imdata.binned.OptoOffTrials = find(sData.behavior.opto.OptoOffTrials>0);
sData.imdata.binned.OptoOnTrials = find(sData.behavior.opto.OptoOnTrials>0); % already corrected for failed trials (discarded)
sData.imdata.binned.AfterOptoTrials = find(sData.behavior.opto.AfterOptoTrials>0);

sData.imdata.binned.ROIsMeanAct_OptoOffTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);
sData.imdata.binned.ROIsMeanAct_OptoOnTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);
sData.imdata.binned.ROIsMeanAct_OptoAfterTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);


%%% preparing data for plots
Xaxis = BinSize:BinSize:BinSize*BinNu;

colors = jet(64);
colors(1,1:3) = 0.5; % set grey when there is no signal

for roi = 1:1:sData.imdata.nROIs
    
    data = sData.imdata.binned.RoidFF{1,roi};
        
    OptoOffdFF = NaN(nTrials,BinNu);
    OptoOndFF = NaN(nTrials,BinNu);
    AfterOptodFF = NaN(nTrials,BinNu);

    OptoOffdFF(sData.imdata.binned.OptoOffTrials,:) = data(sData.imdata.binned.OptoOffTrials,:);
    OptoOndFF(sData.imdata.binned.OptoOnTrials,:) = data(sData.imdata.binned.OptoOnTrials,:);
    AfterOptodFF(sData.imdata.binned.AfterOptoTrials,:) = data(sData.imdata.binned.AfterOptoTrials,:);
          
    sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,:) = nanmean(OptoOffdFF);
    sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,:) = nanmean(OptoOndFF);
    sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,:) = nanmean(AfterOptodFF);

    optoStimStart = sData.behavior.opto.optoStimStart;
    optoStimEnd = sData.behavior.opto.optoStimEnd;
%end
    % for illustration purposes set negative values to zero
    OptoOffdFF(OptoOffdFF<0) = 0;
    OptoOndFF(OptoOndFF<0) = 0;
    AfterOptodFF(AfterOptodFF<0) = 0;

    Cmax = max(max(data));
    if isnan(Cmax)
        continue
    end
    Min = -Cmax/10;
    label = '\DeltaF/F';
     
    %max1 = max(sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu));
    %max2 = max(sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,1:BinNu));
    %max3 = max(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu));
    %Ymax = 1.1*(max([max1,max2,max3])+0.001); 
    YmaxSingleTrials = max(max(sData.imdata.binned.RoidFF{1,roi}));
    Ymax = YmaxSingleTrials;
    Ymax(Ymax==0) = 0.0001;

    figtitle = strcat(sData.sessionInfo.fileID,'-',label,sprintf('-ROI#%d',roi)); 
 
    %%% PLOT
    %optoStimStart = sData.behavior.opto.optoStimStart;
    %optoStimEnd = sData.behavior.opto.optoStimEnd;
    
    figure('Color','white','visible',FigVisible,'Position',[50 50 1000 700]); 
    % Opto-off heatplot
    subplot(3,3,[1,4]);%'Position',[0.1 0.4 0.2 0.4]); %subplot(2,5, i, [l, b, w, h])
    imagesc(Xaxis,1:nTrials,OptoOffdFF); hold on;
    plot(Xaxis,(-(sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)*nTrials)+nTrials),'Color',[1 1 1],'LineWidth',1.5) 
    %line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    colormap(colors); caxis([Min Cmax]); %Cmin/10 for dff  % -Cmax/10 for spikerate
    %xlabel('Position on Wheel (cm)','FontSize',12); 
    ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12; %yticks([1,10,20,30,40]);  ax.FontSize = 12;
    ylabel('Trials','FontSize',12);
    title('Opto-off','FontSize',12);

    % Opto-on heatplot
    subplot(3,3,[2, 5]); 
    imagesc(Xaxis,1:nTrials,OptoOndFF); hold on;
    plot(Xaxis,(-(sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,1:BinNu)*nTrials)+nTrials),'Color',[1 1 1],'LineWidth',1.5) 
    line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    colormap(colors); caxis([Min Cmax]); 
    %xlabel('Position on Wheel (cm)'); 
    ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12;
    ylabel('Trials','FontSize',12);
    title('Opto-on','FontSize',12);
    
    % After-Opto heatplot
    subplot(3,3,[3,6]); 
    imagesc(Xaxis,1:nTrials,AfterOptodFF); hold on;
    plot(Xaxis,(-(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu)*nTrials)+nTrials),'Color',[1 1 1],'LineWidth',1.5) 
    %line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    colormap(colors); caxis([Min Cmax]);
    %xlabel('Position on Wheel (cm)'); 
    ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12;
    ylabel('Trials','FontSize',12);
    title('After-Opto','FontSize',12);
%    sgtitle(figtitle); %suptitle(figtitle); for MAtlab 2018

    % Opto-off mean act
    subplot(3,3,7); 
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    for i = 1:1:nTrials %size(sData.behavior.opto.OptoOffTrialsIndices,1)
        plot(Xaxis,OptoOffdFF(i,1:BinNu),'Color',[0.5 0.5 0.5]); hold on;
    end
    plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu),'Color',[0 0 0],'LineWidth',1) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    %line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12;
    ylabel('Position tuning','FontSize',12);
    %title('Opto-Off');

    % Opto-on mean act
    subplot(3,3,8); 
    %rectangle('Position',[optoStimStart Ymax/100 optoStimEnd-optoStimStart Ymax*0.95],'FaceColor',[1 0.9 0.9],'EdgeColor','none','Alpha',0.5);
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
      for i = 1:1:nTrials %size(sData.behavior.opto.OptoOnTrialsIndices,1)
         plot(Xaxis,OptoOndFF(i,1:BinNu),'Color',[0.5 0.5 0.5]); hold on;
      end
    plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,1:BinNu),'Color',[0 0 0],'LineWidth',1)  %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOnTrials(roi,1:BinNu)
    line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]);
    xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12;
    ylabel('Position tuning','FontSize',12);
    text(BinSize*2,Ymax*0.9,'638 nm laser is on','Color','red','FontSize',12);
    %title('Opto-On');

    % after-Opto mean act
    subplot(3,3,9);
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 nTrials],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
    for i = 1:1:nTrials 
        plot(Xaxis,AfterOptodFF(i,1:BinNu),'Color',[0.5 0.5 0.5]); hold on;
    end
    plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu),'Color',[0 0 0],'LineWidth',1)   %Xaxis,CADATA.Opto.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu)
    %line([optoStimStart optoStimStart],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 nTrials],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]);
    xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 12;
    ylabel('Position tuning','FontSize',12);
    %title('After-Opto');
    FileName = strcat(sprintf('ROI-%d',roi),'-','-opto'); 
    savefig(fullfile(SavePath,FileName));
    saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));
 
    if rem(roi,20) == 0
           close all;
    end
end 
close all;


Ymax = 1.2*max([nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1),nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1),nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1)]);

%mean activity of all ROIs Opto-off-on-after
figure('Color','white','Position',[100 100 800 250]);
subplot(1,3,1);
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1),'LineWidth',1.5) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
%line([optoStimStart optoStimStart],[0 Ymax],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1); hold on; 
axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)','FontSize',10); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 10;
ylabel('Position tuning','FontSize',10);
title('Opto-Off','FontSize',12);

subplot(1,3,2);
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1),'LineWidth',1.5) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
line([optoStimStart optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1.5); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)','FontSize',10); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 10;
%ylabel('Position tuning','FontSize',12);
text(BinSize*3,Ymax*0.9,'638 nm laser is on','Color','red','FontSize',10);
title('Opto-On','FontSize',12);

subplot(1,3,3);
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); hold on; 
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 Ymax],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1),'LineWidth',1.5) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
%line([optoStimStart optoStimStart],[0 Ymax],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color',[0.5 0.5 0.5],'LineStyle','-','LineWidth',1); hold on; 
axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)','FontSize',10); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]); ax.FontSize = 10;
%ylabel('Position tuning','FontSize',12);
title('After-Opto','FontSize',12);

FileName = 'AllROIMeanAct-opto'; 
savefig(fullfile(SavePath,FileName));
saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));

%mean activity of all ROIs Opto-off-on-after in one plot
figure('Color','white','Position',[100 100 400 300]);
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1),'LineStyle','-','LineWidth',1.5) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
hold on
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1),'LineStyle','-','LineWidth',1.5)
hold on
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1),'LineStyle','-','LineWidth',1)
hold on
line([optoStimStart optoStimStart],[0 Ymax*1.2],'Color','red','LineStyle','-','LineWidth',1.5); hold on; line([optoStimEnd optoStimEnd],[0 Ymax*1.2],'Color','red','LineStyle','-','LineWidth',1); 
line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 Ymax*1.2],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 Ymax*1.2],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5);  
line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 Ymax*1.2],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 Ymax*1.2],'Color',[0.8 0.8 0.8],'LineStyle','--','LineWidth',1.5); 
legend('Opto-off', 'Opto-on', 'After-opto','Location','Northeast');
axis([0 160 0 Ymax*1.2]); % ceil(Ymax)
xlabel('Position on Wheel (cm)','FontSize',12); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning','FontSize',12);
text(BinSize*3,Ymax*1.1,'638 nm laser is on','Color','red','FontSize',12);
title('Position tuning of all cells','FontSize',14);
FileName = 'AllROIMeanAct-opto2'; 
savefig(fullfile(SavePath,FileName));
saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end