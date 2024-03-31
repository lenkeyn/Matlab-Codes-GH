function sData = PCFirstSecondHalf(sData,LapsTested)

% LapsTested = 50;

sData.imdata.MaoPC.FirstSecondHalf = struct;
%savePath = 'C:\MATLAB\SAVE';
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'PlaceCell\MeanActivitySortedFirstSecondHalf');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\PlaceCell\MeanActivitySortedFirstSecondHalf');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
Yaxis = 1:1:length(sData.imdata.MaoPC.AllPlaceCells); 

% set which laps belong to the first half and second half of the session 
LastLap = sData.behavior.wheelLapImaging-1; 
%LapsTested = 50; 
FirstStartLap = 1;
FirstEndLap = FirstStartLap + LapsTested - 1;
SecondStartLap = LastLap - LapsTested + 1;
SecondEndLap = LastLap;

sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsFirst = sData.imdata.binned.LightOffTrials(find(sData.imdata.binned.LightOffTrials >= FirstStartLap,1):find(sData.imdata.binned.LightOffTrials <= FirstEndLap,1,'last'));
sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsSecond = sData.imdata.binned.LightOffTrials(find(sData.imdata.binned.LightOffTrials >= SecondStartLap,1):find(sData.imdata.binned.LightOffTrials <= SecondEndLap,1,'last'));
sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsFirst = sData.imdata.binned.LightOnTrials(find(sData.imdata.binned.LightOnTrials >= FirstStartLap,1):find(sData.imdata.binned.LightOnTrials <= FirstEndLap,1,'last'));
sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsSecond = sData.imdata.binned.LightOnTrials(find(sData.imdata.binned.LightOnTrials >= SecondStartLap,1):find(sData.imdata.binned.LightOnTrials <= SecondEndLap,1,'last'));

sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffFirst = NaN(length(sData.imdata.MaoPC.AllPlaceCells),sData.behavior.meta.nBins);
sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffSecond = NaN(length(sData.imdata.MaoPC.AllPlaceCells),sData.behavior.meta.nBins);

for i = 1:1:length(sData.imdata.MaoPC.AllPlaceCells) % collect each place roi data light-on  light-off trials only in first or second part
    j = sData.imdata.MaoPC.AllPlaceCells(i); % j is the original ROI name
    sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffFirst(i,:) = nanmean(sData.imdata.binned.RoidFF{1,j}(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsFirst,:));  
    sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnFirst(i,:) = nanmean(sData.imdata.binned.RoidFF{1,j}(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsFirst,:));  
    sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffSecond(i,:) = nanmean(sData.imdata.binned.RoidFF{1,j}(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsSecond,:));  
    sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnSecond(i,:) = nanmean(sData.imdata.binned.RoidFF{1,j}(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsSecond,:));  
end

%%% calculating max value among laser on and laser off trials
Max = NaN(length(sData.imdata.MaoPC.AllPlaceCells),4);
Max(:,1) = max(sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffFirst,[],2);
Max(:,2) = max(sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnFirst,[],2);
Max(:,3) = max(sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffFirst,[],2);
Max(:,4) = max(sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnFirst,[],2);
MaxAct = max(Max,[],2);


%%% PLOT
% Laser-Off, first half
HeatMapData = sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffFirst; 
NormHeatMapData = HeatMapData./MaxAct;
SortedData1 = NormHeatMapData(sData.imdata.MaoPC.SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:length(sData.imdata.MaoPC.AllPlaceCells); yticks = linspace(1, length(sData.imdata.MaoPC.AllPlaceCells), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat('Place cell activity - Laser-off - lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap))); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-dff-Sorted-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap)); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Laser-On, first half
HeatMapData = sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnFirst; 
NormHeatMapData = HeatMapData./MaxAct;
SortedData1 = NormHeatMapData(sData.imdata.MaoPC.SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:length(sData.imdata.MaoPC.AllPlaceCells); yticks = linspace(1, length(sData.imdata.MaoPC.AllPlaceCells), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat('Place cell activity - Laser-on - lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap))); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-dff-SortedMeanAct-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap)); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%%% PLOT
% Laser-Off, second half
HeatMapData = sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOffSecond; 
NormHeatMapData = HeatMapData./MaxAct;
SortedData1 = NormHeatMapData(sData.imdata.MaoPC.SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:length(sData.imdata.MaoPC.AllPlaceCells); yticks = linspace(1, length(sData.imdata.MaoPC.AllPlaceCells), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat('Place cell activity - Laser-off - lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap))); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOff-dff-Sorted-lap-',num2str(SecondStartLap),'-',num2str(SecondEndLap)); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Laser-On, second half
HeatMapData = sData.imdata.MaoPC.FirstSecondHalf.MeanROIActLightOnSecond; 
NormHeatMapData = HeatMapData./MaxAct;
SortedData1 = NormHeatMapData(sData.imdata.MaoPC.SortingOrder(:,2),:); % new matrix containing sorted data.
%PLOT FIGURE
figure('Color','white'); 
imagesc(Xaxis,Yaxis,SortedData1) %(1:number of bins;1:number of trials)
c = colorbar; colormap(jet); caxis([0 1]);
c.Label.String = 'dFF'; c.Label.FontSize = 11; c.TickDirection = 'out'; 
line([sData.behavior.opto.optoStimStart sData.behavior.opto.optoStimStart],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2); hold on;    
line([sData.behavior.opto.optoStimEnd sData.behavior.opto.optoStimEnd],[0 length(sData.imdata.MaoPC.AllPlaceCells)],'Color','white','LineStyle','-','LineWidth',2);    
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
ylabel('ROIs'); yticklabels = 0:10:length(sData.imdata.MaoPC.AllPlaceCells); yticks = linspace(1, length(sData.imdata.MaoPC.AllPlaceCells), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(strcat('Place cell activity - Laser-on - lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap))); 
FileName = strcat(sData.sessionInfo.fileID,'-LaserOn-dff-SortedMeanAct-lap-',num2str(SecondStartLap),'-',num2str(SecondEndLap)); 
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

%{
%%% PLOT mean Licks
sData.behavior.FirstSecondHalf = struct;
sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffFirst = sData.behavior.binning.lickPerCmBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsFirst,:);
sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnFirst = sData.behavior.binning.lickPerCmBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsFirst,:);
sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffSecond = sData.behavior.binning.lickPerCmBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsSecond,:);
sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnSecond = sData.behavior.binning.lickPerCmBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsSecond,:);

% mean Lick First half
RewardZoneLength = 6;
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*(sData.behavior.meta.nBins+10);
Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffFirst)));
rectangle('Position',[sData.behavior.opto.optoStimStart Ymax/500 sData.behavior.opto.optoStimEnd-sData.behavior.opto.optoStimStart Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffFirst)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnFirst)); hold on;
line([157 157],[0 Ymax],'Color','black','LineStyle','--'); hold on;
line([157+RewardZoneLength 157+RewardZoneLength],[0 Ymax],'Color','black','LineStyle','--');
xlabel('Position on wheel (cm)');
ylabel('Licks/cm');
legend(strcat('Laser-Off lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Laser-On lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),'Location','North');
FileName = strcat(sData.sessionInfo.fileID,'LickFirstHalf-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% mean Lick Second half
RewardZoneLength = 6;
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*(sData.behavior.meta.nBins+10);
Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffSecond)));
rectangle('Position',[sData.behavior.opto.optoStimStart Ymax/500 sData.behavior.opto.optoStimEnd-sData.behavior.opto.optoStimStart Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffSecond)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnSecond)); hold on;
line([157 157],[0 Ymax],'Color','black','LineStyle','--'); hold on;
line([157+RewardZoneLength 157+RewardZoneLength],[0 Ymax],'Color','black','LineStyle','--');
xlabel('Position on wheel (cm)');
ylabel('Licks/cm');
legend(strcat('Laser-Off lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),strcat('Laser-On lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),'Location','North');
FileName = strcat(sData.sessionInfo.fileID,'LickSecondHalf-lap-',num2str(SecondStartLap),'-',num2str(SecondEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% mean Lick First-Second half laser on-off
RewardZoneLength = 6;
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*(sData.behavior.meta.nBins+10);
Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffFirst)));
rectangle('Position',[sData.behavior.opto.optoStimStart Ymax/500 sData.behavior.opto.optoStimEnd-sData.behavior.opto.optoStimStart Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffFirst)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnFirst)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOffSecond)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.LickPerCmBinned.LightOnSecond)); hold on;
line([157 157],[0 Ymax],'Color','black','LineStyle','--'); hold on;
line([157+RewardZoneLength 157+RewardZoneLength],[0 Ymax],'Color','black','LineStyle','--');
xlabel('Position on wheel (cm)');
ylabel('Licks/cm');
legend(strcat('Laser-Off lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Laser-On lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Laser-Off lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),strcat('Laser-On lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),'Location','North');
FileName = strcat(sData.sessionInfo.fileID,'LickFirstSecondHalf-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap),'-',num2str(SecondStartLap),'-',num2str(SecondEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));


%%% PLOT mean Velo
sData.behavior.FirstSecondHalf.VeloBinned.LightOffFirst = sData.behavior.binning.veloBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsFirst,:);
sData.behavior.FirstSecondHalf.VeloBinned.LightOnFirst = sData.behavior.binning.veloBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsFirst,:);
sData.behavior.FirstSecondHalf.VeloBinned.LightOffSecond = sData.behavior.binning.veloBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOffTrialsSecond,:);
sData.behavior.FirstSecondHalf.VeloBinned.LightOnSecond = sData.behavior.binning.veloBinned(sData.imdata.MaoPC.FirstSecondHalf.LightOnTrialsSecond,:);

% mean Velo First half
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOffFirst)));
rectangle('Position',[sData.behavior.opto.optoStimStart Ymax/500 sData.behavior.opto.optoStimEnd-sData.behavior.opto.optoStimStart Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOffFirst)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOnFirst)); hold on;
xlabel('Position on wheel (cm)');
ylabel('Speed (cm/s)');
legend(strcat('Laser-Off lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),strcat('Laser-On lap:',num2str(FirstStartLap),'-',num2str(FirstEndLap)),'Location','South');
FileName = strcat(sData.sessionInfo.fileID,'SpeedFirstHalf-lap-',num2str(FirstStartLap),'-',num2str(FirstEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% mean Velo Second half
figure('Color','white'); 
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
Ymax = ceil(max(mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOffSecond)));
rectangle('Position',[sData.behavior.opto.optoStimStart Ymax/500 sData.behavior.opto.optoStimEnd-sData.behavior.opto.optoStimStart Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOffSecond)); hold on;
plot(Xaxis,mean(sData.behavior.FirstSecondHalf.VeloBinned.LightOnSecond)); hold on;
xlabel('Position on wheel (cm)');
ylabel('Speed (cm/s)');
legend(strcat('Laser-Off lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),strcat('Laser-On lap:',num2str(SecondStartLap),'-',num2str(SecondEndLap)),'Location','South');
FileName = strcat(sData.sessionInfo.fileID,'SpeedSecondHalf-lap-',num2str(SecondStartLap),'-',num2str(SecondEndLap));
savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
%}

% Save file to same path where other files can be found 
%save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end


