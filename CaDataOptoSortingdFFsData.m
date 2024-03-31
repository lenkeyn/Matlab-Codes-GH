function sData  = CaDataOptoSortingdFFsData(sData,datatype,IsOpto) 

% type: which data type to use, it can be 0 = RoidFF or 1 = RoidFF-SR (slow transients removed), 2= RoidFF-SR-LP (slow transients removed and transients are low pass filtered)
% IsOpto  = 1 : it was an optically stimulated session, IsOpto = 0 : sessions without stimulation

mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiAct-Opto');
SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct-Opto');
%SavePath = sData.sessionInfo.savePath;
BinSize = sData.behavior.meta.binSize;
BinNu = sData.behavior.meta.nBins;
TRNu = sData.behavior.wheelLapImaging;

%Sorting binned trials for plotting
sData.imdata.binned.OptoOffTrials = find(sData.behavior.opto.OptoOffTrials>0);
sData.imdata.binned.OptoOnTrials = find(sData.behavior.opto.OptoOnTrials>0);
sData.imdata.binned.AfterOptoTrials = find(sData.behavior.opto.AfterOptoTrials>0);

sData.imdata.binned.ROIsMeanAct_OptoOffTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);
sData.imdata.binned.ROIsMeanAct_OptoOnTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);
sData.imdata.binned.ROIsMeanAct_OptoAfterTrials = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);


%%% preparing data for plots
Xaxis = BinSize:BinSize:BinSize*BinNu;

colors = jet(64);
colors(1,1:3) = 0.5; % set grey when there is no signal

 %%% calculate opto-stim
    %{
    optoStimMax = max(mean(sData.behavior.opto.OptoStimOnMatrix),[],2);
    meanOptoStim = zeros(1,sData.behavior.meta.nBins);
    meanOptoStim((mean(sData.behavior.opto.OptoStimOnMatrix)) > optoStimMax/4) = 1;
    optoStimStart = find(diff(meanOptoStim) == 1)*sData.behavior.meta.binSize;
    optoStimEnd = find(diff(meanOptoStim) == -1)*sData.behavior.meta.binSize;
    if  isempty(optoStimStart) 
        optoStimStart = 5;
    end
    if  isempty(optoStimEnd) 
        optoStimEnd = 157;
    elseif optoStimEnd < optoStimStart
        optoStimEnd = 157;
    end
    sData.behavior.opto.optoStimStart = optoStimStart;
    sData.behavior.opto.optoStimEnd = optoStimEnd;
    %}



for roi = 1:1:sData.imdata.nROIs
    
    if datatype == 0
        data = sData.imdata.binned.RoidFF{1,roi};
    elseif datatype == 1
        data = sData.imdata.binned.RoidFF_SR{1,roi};
    elseif datatype == 2
        data = sData.imdata.binned.RoidFF_SR_LP{1,roi};    
    end

    OptoOffdFF = NaN(sData.behavior.wheelLapImaging,BinNu);
    OptoOndFF = NaN(sData.behavior.wheelLapImaging,BinNu);
    AfterOptodFF = NaN(sData.behavior.wheelLapImaging,BinNu);

    OptoOffdFF((sData.behavior.opto.OptoOffTrials>0),:) = data((sData.behavior.opto.OptoOffTrials>0),:);
    OptoOndFF((sData.behavior.opto.OptoOnTrials>0),:) = data((sData.behavior.opto.OptoOnTrials>0),:);
    AfterOptodFF((sData.behavior.opto.AfterOptoTrials>0),:) = data((sData.behavior.opto.AfterOptoTrials>0),:);
    
      
    sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,:) = nanmean(OptoOffdFF);
    sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,:) = nanmean(OptoOndFF);
    sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,:) = nanmean(AfterOptodFF);


    Cmax = max(max(data));
    Cmin = min(min(data));
    text = 'dFF';
     
    max1 = max(sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu));
    max2 = max(sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,1:BinNu));
    max3 = max(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu));
    Ymax = 1.1*(max([max1,max2,max3])+0.001);

    figtitle = strcat(sData.sessionInfo.fileID,'-',text,sprintf('-ROI#%d',roi)); 

    if IsOpto == 1 
        
        %%% PLOT
        figure('Color','white','Position',[100 1000 1000 700]); 
        % Opto-off heatplot
        subplot(3,3,[1,4]);%'Position',[0.1 0.4 0.2 0.4]); %subplot(2,5, i, [l, b, w, h])
        imagesc(Xaxis,1:TRNu,OptoOffdFF); %Xaxis,1:TRNu,OptoOffdFF
        line([optoStimStart optoStimStart],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; 
        %c = colorbar; c.TickDirection = 'out';
        colormap(colors); caxis([Cmin Cmax]); %Cmin/10 for dff  % -Cmax/10 for spikerate
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials');
        title('Laser-off');

        % Opto-on heatplot
        %figure('Color','white','Position',[100 100 1000 600]); 
        subplot(3,3,[2, 5]); %,'Position',[0.4 0.4 0.2 0.4]   %subplot(3,3,[2,5]);
        imagesc(Xaxis,1:TRNu,OptoOndFF); hold on; %Xaxis,1:TRNu,OptoOndFF
        line([optoStimStart optoStimStart],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; 
        %c = colorbar;  c.TickDirection = 'out'; 
        colormap(colors); caxis([Cmin Cmax]); 
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials');
        title('Laser-on');

        % Opto-on heatplot
        %figure('Color','white','Position',[100 100 1000 600]); 
        subplot(3,3,[3, 6]); %,'Position',[0.7 0.4 0.2 0.4]);  %subplot(3,3,[3,6]);
        imagesc(Xaxis,1:TRNu,AfterOptodFF);  %Xaxis,1:TRNu,AfterOptodFF
        line([optoStimStart optoStimStart],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','white','LineStyle','-','LineWidth',1); hold on; 
        %c = colorbar; c.TickDirection = 'out'; 
        colormap(colors); caxis([Cmin Cmax]);
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials');
        title('After-Laser');
        suptitle(figtitle);

        % Opto-off mean act
        subplot(3,3,7); 
        plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
        line([optoStimStart optoStimStart],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        axis([0 160 0 Ymax]); % ceil(Ymax)
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning');
        title('Laser-Off');

        % Opto-on mean act
        subplot(3,3,8); 
        plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoOnTrials(roi,1:BinNu))  %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOnTrials(roi,1:BinNu)
        line([optoStimStart optoStimStart],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        axis([0 160 0 Ymax]);
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning');
        title('Laser-On');

        % after-Opto mean act
        subplot(3,3,9); 
        plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu))   %Xaxis,CADATA.Opto.ROIsMeanAct_OptoAfterTrials(roi,1:BinNu)
        line([optoStimStart optoStimStart],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 TRNu],'Color','red','LineStyle','-','LineWidth',1); hold on; 
        axis([0 160 0 Ymax]);
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning');
        title('After-Laser');
        FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-opto'); 
        savefig(fullfile(SavePath,FileName));
        saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));
        
        
    
    elseif IsOpto == 0
        
        %%% PLOT
        figure('Color','white','Position',[100 100 500 1000]); 
        % Opto-off heatplot
        subplot(3,1,[1,2]); %,'Position',[0.1 0.4 0.2 0.4]); %subplot(2,5, i, [l, b, w, h])
        imagesc(Xaxis,1:TRNu,OptoOffdFF);
        %c = colorbar; c.TickDirection = 'out';
        colormap(colors); caxis([Cmin Cmax]); %Cmin/10 for dff  % -Cmax/10 for spikerate
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Trials');
        suptitle(figtitle);
        
        % Opto-off mean act
        subplot(3,1,3); 
        plot(Xaxis,sData.imdata.binned.ROIsMeanAct_OptoOffTrials(roi,1:BinNu))
        axis([0 160 0 Ymax]); % ceil(Ymax)
        xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
        ylabel('Position tuning');
        FileName = strcat(sprintf('ROI-%d',roi),'-',text,'-opto'); 
        savefig(fullfile(SavePath,FileName));
        saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));
              
        
    end

end

close all;

if IsOpto == 1
    Ymax = 1.1*max([nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1),nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1),nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1)]);
    
    %mean activity of all ROIs Opto-off-on-after
    figure('Color','white','Position',[100 100 1500 300]);
    subplot(1,3,1);
    plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1)) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    line([optoStimStart optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('Laser-Off');

    subplot(1,3,2);
    plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1)) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    line([optoStimStart optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('Laser-On');

    subplot(1,3,3);
    plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1)) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
    line([optoStimStart optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
    axis([0 160 0 Ymax]); % ceil(Ymax)
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning');
    title('After-Laser');

    FileName = 'AllROIMeanAct-opto'; 
    savefig(fullfile(SavePath,FileName));
    saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));
end

%mean activity of all ROIs Opto-off-on-after in one plot
figure('Color','white','Position',[100 100 500 400]);
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOffTrials,1),'LineStyle','-','LineWidth',1) %Xaxis,CADATA.Opto.ROIsMeanAct_OptoOffTrials(roi,1:BinNu)
hold on
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoOnTrials,1),'LineStyle','-','LineWidth',1)
hold on
plot(Xaxis,nanmean(sData.imdata.binned.ROIsMeanAct_OptoAfterTrials,1),'LineStyle','-','LineWidth',1)
hold on
line([optoStimStart optoStimStart],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; line([optoStimEnd optoStimEnd],[0 Ymax],'Color','red','LineStyle','-','LineWidth',1); hold on; 
axis([0 160 0 Ymax]); % ceil(Ymax)
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning');
legend('Laser-off', 'Laser-on', 'After-laser','Location','South');
title('Position tuning of all cells');
FileName = 'AllROIMeanAct-opto2'; 
savefig(fullfile(SavePath,FileName));
saveas(gcf,(fullfile(SavePath,[FileName '.jpg'])));


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end