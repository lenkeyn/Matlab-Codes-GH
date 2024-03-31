function PlaceCellComp = PosTuningComp(sData1,sData2,type)
% comparing pos tuning curves in same ROIs in different sessions

% type = 0; type 0 is dff, type is deconvolved signal
% sData1;
% sData2;
SessionName1 = 'S1';
SessionName2 = 'S2';

nBins = sData1.behavior.meta.nBins;
BinSize = sData1.behavior.meta.binSize;
nROIs = sData1.imdata.nROIs;

if type == 0 
    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PosTuning-comparison-dff-PCEither');
    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PosTuning-comparison-dff-PCBoth');
    MaoPC1 = sData1.imdata.MaoPC_dff;
    MaoPC2 = sData2.imdata.MaoPC_dff;
elseif type == 1
    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PosTuning-comparison-deconv-PCEither');
    mkdir(strcat(sData1.sessionInfo.savePath,'\Imaging'),'PosTuning-comparison-deconv-PCBoth');
    MaoPC1 = sData1.imdata.MaoPC_deconv;
    MaoPC2 = sData2.imdata.MaoPC_deconv;
end

PlaceCellComp = struct;
PlaceCellComp.EitherPlaceCells = struct;
PlaceCellComp.BothPlaceCells = struct;
% Collect place cells which is a place cell in either S1 or S2
EitherPlaceCellsPre = NaN(nROIs,1);
BothPlaceCellsPre = NaN(nROIs,1);
for i=1:1:nROIs
    if sum(any(MaoPC1.PlaceCells(:,1)==i)) || sum(any(MaoPC2.PlaceCells(:,1)==i)) %sData.imdata.MaoPC.LightOff.PlaceCells(:,1)==i    sData1.imdata.MaoPC_dff.PlaceCells
       EitherPlaceCellsPre(i,1) = i; 
    end
    if sum(any(MaoPC1.PlaceCells(:,1)==i)) && sum(any(MaoPC2.PlaceCells(:,1)==i)) %sData.imdata.MaoPC.LightOff.PlaceCells(:,1)==i    sData1.imdata.MaoPC_dff.PlaceCells
       BothPlaceCellsPre(i,1) = i; 
    end
end
EitherPlaceCellsROIs = zeros(length(EitherPlaceCellsPre(~isnan(EitherPlaceCellsPre))),3);
EitherPlaceCellsROIs(:,1) = EitherPlaceCellsPre(~isnan(EitherPlaceCellsPre));
BothPlaceCellsROIs = zeros(length(BothPlaceCellsPre(~isnan(BothPlaceCellsPre))),3);
BothPlaceCellsROIs(:,1) = BothPlaceCellsPre(~isnan(BothPlaceCellsPre));
PlaceCellComp.BothPlaceCells.ROIs(:,1) = BothPlaceCellsROIs(:,1);
nBothPC = length(BothPlaceCellsROIs);

% write if the place cell was a place cell in either S1 or S2
for i=1:1:length(EitherPlaceCellsROIs)
   if MaoPC1.Criteria123Passed(EitherPlaceCellsROIs(i,1))== EitherPlaceCellsROIs(i,1)   
       EitherPlaceCellsROIs(i,2) = 1; 
   end  
   if MaoPC2.Criteria123Passed(EitherPlaceCellsROIs(i,1))== EitherPlaceCellsROIs(i,1)
       EitherPlaceCellsROIs(i,3) = 1; 
   end  
end
PlaceCellComp.EitherPlaceCells.ROIs = EitherPlaceCellsROIs(:,1);
nEitherPC = length(EitherPlaceCellsROIs);

% what is in the columns
columns(1,1) = {'ROIid'}; columns(1,2) = {'S1'}; columns(1,3) = {'S2'}; columns(1,4) = {'difference'};
PlaceCellComp.columns = columns;

% calculate peak of pos tuning curve shift cell shift in all place cells (either in S1 or S2), search for the peak (which bin)
PlaceCellComp.EitherPlaceCells.PosTuningPeak = zeros(nEitherPC,4);
PlaceCellComp.EitherPlaceCells.PosTuningPeak(:,1) = EitherPlaceCellsROIs(:,1); 
for i=1:1:nEitherPC
    PlaceCellComp.EitherPlaceCells.PosTuningPeak(i,2) = find(MaoPC1.PosTuningNorm(EitherPlaceCellsROIs(i,1),:)==1);
    PlaceCellComp.EitherPlaceCells.PosTuningPeak(i,3) = find(MaoPC2.PosTuningNorm(EitherPlaceCellsROIs(i,1),:)==1);
    PlaceCellComp.EitherPlaceCells.PosTuningPeak(i,4) = PlaceCellComp.EitherPlaceCells.PosTuningPeak(i,3) - PlaceCellComp.EitherPlaceCells.PosTuningPeak(i,2);
end
% calculate place cell shift if the place cell remained, search for the peak (which bin)
PlaceCellComp.BothPlaceCells.PosTuningPeak = zeros(nBothPC,4);
PlaceCellComp.BothPlaceCells.PosTuningPeak(:,1) = BothPlaceCellsROIs(:,1); 
for i=1:1:nBothPC
    PlaceCellComp.BothPlaceCells.PosTuningPeak(i,2) = find(MaoPC1.PosTuningNorm(BothPlaceCellsROIs(i,1),:)==1);
    PlaceCellComp.BothPlaceCells.PosTuningPeak(i,3) = find(MaoPC2.PosTuningNorm(BothPlaceCellsROIs(i,1),:)==1);
    PlaceCellComp.BothPlaceCells.PosTuningPeak(i,4) = PlaceCellComp.BothPlaceCells.PosTuningPeak(i,3) - PlaceCellComp.BothPlaceCells.PosTuningPeak(i,2);
end

% calculate place field widening or narrowing
PlaceCellComp.EitherPlaceCells.PlaceFieldWidth = zeros(nEitherPC,4);
PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(:,1) = EitherPlaceCellsROIs(:,1); 
PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(:,2) = 0;
PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(:,3) = 0;
for i=1:1:nEitherPC
    Temp1 = MaoPC1.PlaceFieldBinLength(EitherPlaceCellsROIs(i,1),:);
    Temp2 = MaoPC2.PlaceFieldBinLength(EitherPlaceCellsROIs(i,1),:);
    if ~isnan(Temp1(~isnan(Temp1)))
        PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(i,2) = max(Temp1(~isnan(Temp1)));
    end
    if ~isnan(Temp2(~isnan(Temp2)))
        PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(i,3) = max(Temp2(~isnan(Temp2)));
    end
    PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(i,4) = PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(i,3) - PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(i,2);
end
% calculate place cell shift if the place cell remained, search for the peak (which bin)
PlaceCellComp.BothPlaceCells.PlaceFieldWidth = zeros(nBothPC,4);
PlaceCellComp.BothPlaceCells.PlaceFieldWidth(:,1) = BothPlaceCellsROIs(:,1); 
PlaceCellComp.BothPlaceCells.PlaceFieldWidth(:,2) = 0;
PlaceCellComp.BothPlaceCells.PlaceFieldWidth(:,3) = 0;
for i=1:1:nBothPC
    Temp1 = MaoPC1.PlaceFieldBinLength(BothPlaceCellsROIs(i,1),:);
    Temp2 = MaoPC2.PlaceFieldBinLength(BothPlaceCellsROIs(i,1),:);
    if ~isnan(Temp1(~isnan(Temp1)))
        PlaceCellComp.BothPlaceCells.PlaceFieldWidth(i,2) = max(Temp1(~isnan(Temp1)));
    end
    if ~isnan(Temp2(~isnan(Temp2)))
        PlaceCellComp.BothPlaceCells.PlaceFieldWidth(i,3) = max(Temp2(~isnan(Temp2)));
    end
    PlaceCellComp.BothPlaceCells.PlaceFieldWidth(i,4) = PlaceCellComp.BothPlaceCells.PlaceFieldWidth(i,3) - PlaceCellComp.BothPlaceCells.PlaceFieldWidth(i,2);
end

% calculate place cell gain/loss S1 vs S2
PlaceCellComp.PlaceCellMoreS2vsS1 = MaoPC2.placeROINu - MaoPC1.placeROINu;
PlaceCellComp.PlaceCellMorePercentS2vsS1 = PlaceCellComp.PlaceCellMoreS2vsS1/MaoPC2.placeROINu*100;
% mean shift in pos tuning
PlaceCellComp.MeanPosTuningShiftEither = mean(PlaceCellComp.EitherPlaceCells.PosTuningPeak(:,4));
PlaceCellComp.MeanPosTuningShiftBoth = mean(PlaceCellComp.BothPlaceCells.PosTuningPeak(:,4));
% mean place field width change 
PlaceCellComp.MeanPlaceFieldWidthChangeEither = mean(PlaceCellComp.EitherPlaceCells.PlaceFieldWidth(:,4));
PlaceCellComp.MeanPlaceFieldWidthChangeBoth = mean(PlaceCellComp.BothPlaceCells.PlaceFieldWidth(:,4));


% plot place cells pos tuning curve in both sessions (either place cell in S1 or S2):
if type == 0
    SavePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PosTuning-comparison-dff-PCEither');
elseif type == 1
    SavePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PosTuning-comparison-deconv-PCEither');
end
for j = 1:1:nEitherPC
    i = PlaceCellComp.EitherPlaceCells.ROIs(j);
    figure();
    Xaxis = BinSize:BinSize:BinSize*nBins;
    Ymax = max(max(MaoPC1.PosTuning(i,:)),max(MaoPC2.PosTuning(i,:)))+0.00000000001;
    plot(Xaxis,MaoPC1.PosTuning(i,1:nBins),'LineWidth',2); hold on; % pos tuning curve
    plot(Xaxis,MaoPC2.PosTuning(i,1:nBins),'LineWidth',2); hold on; % pos tuning curve
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData1.sessionInfo.fileID,' ROI #',num2str(i),' Pos. tuning'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    legend(SessionName1,SessionName2);
    fname = strcat(sData1.sessionInfo.fileID,'-roi',num2str(i),'_PosTuning');
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
 end
close all

% plot mean
PlaceCellComp.MeanPosTuningEitherS1 = mean(MaoPC1.PosTuning(PlaceCellComp.EitherPlaceCells.ROIs,:),1);
PlaceCellComp.MeanPosTuningEitherS2 = mean(MaoPC2.PosTuning(PlaceCellComp.EitherPlaceCells.ROIs,:),1);
figure();
Ymax = max(max(PlaceCellComp.MeanPosTuningEitherS1),max(PlaceCellComp.MeanPosTuningEitherS2))+0.0000001;
plot(Xaxis,PlaceCellComp.MeanPosTuningEitherS1,'LineWidth',2); hold on; % pos tuning curve
plot(Xaxis,PlaceCellComp.MeanPosTuningEitherS2,'LineWidth',2); hold on; % pos tuning curve
axis([0 160 0 Ymax]); % ceil(Ymax)
title(strcat(sData1.sessionInfo.fileID,' Pos. tuning MEAN either'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
legend(SessionName1,SessionName2);
fname = strcat(sData1.sessionInfo.fileID,'-roi',num2str(i),'_PosTuning-MEAN(either)');
savefig(fullfile(SavePath,fname));
saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));


% plot place cells pos tuning curve in both sessions (either in S1 or S2):
if type == 0
    SavePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PosTuning-comparison-dff-PCBoth');
elseif type == 1
    SavePath = strcat(sData1.sessionInfo.savePath,'\Imaging\PosTuning-comparison-deconv-PCBoth');
end
for j = 1:1:nBothPC
    i = PlaceCellComp.BothPlaceCells.ROIs(j);
    figure();
    Xaxis = BinSize:BinSize:BinSize*nBins;
    Ymax = max(max(MaoPC1.PosTuning(i,:)),max(MaoPC2.PosTuning(i,:)))+0.00000000001;
    %rectangle('Position',[MaoPC1.PlaceFieldStartBin(i,Column(i))*BinSize Ymax/500 MaoPC1.PlaceFieldBinLength(i,Column(i))*BinSize Ymax],'FaceColor',[1 0.95 1],'EdgeColor','none'); hold on;
    plot(Xaxis,MaoPC1.PosTuning(i,1:nBins),'LineWidth',2); hold on; % pos tuning curve
    plot(Xaxis,MaoPC2.PosTuning(i,1:nBins),'LineWidth',2); hold on; % pos tuning curve
    %line([0 BinSize*nBins],[Treshold(i) Treshold(i)],'Color','red','LineWidth',1); hold on; 
    line([PlaceCellComp.BothPlaceCells.PosTuningPeak(j,2)*BinSize PlaceCellComp.BothPlaceCells.PosTuningPeak(j,2)*BinSize],[0 Ymax],'Color',[0 0.45 0.74],'LineStyle','--','LineWidth',2); hold on; 
    line([PlaceCellComp.BothPlaceCells.PosTuningPeak(j,3)*BinSize PlaceCellComp.BothPlaceCells.PosTuningPeak(j,3)*BinSize],[0 Ymax],'Color','red','LineStyle','--','LineWidth',2); hold on; 
    line([MaoPC1.PotPlaceFieldPos(i,1)*BinSize MaoPC1.PotPlaceFieldPos(i,1)*BinSize],[0 Ymax],'Color',[0 0.45 0.74],'LineStyle',':','LineWidth',0.5); hold on; line([(MaoPC1.PotPlaceFieldPos(i,1)+MaoPC1.PotPlaceFieldLength(i,1))*BinSize (MaoPC1.PotPlaceFieldPos(i,1)+MaoPC1.PotPlaceFieldLength(i,1))*BinSize],[0 Ymax],'Color',[0 0.45 0.74],'LineStyle',':','LineWidth',0.5); hold on;
    line([MaoPC2.PotPlaceFieldPos(i,1)*BinSize MaoPC2.PotPlaceFieldPos(i,1)*BinSize],[0 Ymax],'Color','red','LineStyle',':','LineWidth',0.5); hold on; line([(MaoPC2.PotPlaceFieldPos(i,1)+MaoPC2.PotPlaceFieldLength(i,1))*BinSize (MaoPC2.PotPlaceFieldPos(i,1)+MaoPC2.PotPlaceFieldLength(i,1))*BinSize],[0 Ymax],'Color','red','LineStyle',':','LineWidth',0.5); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData1.sessionInfo.fileID,' ROI #',num2str(i),' Pos. tuning'));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    legend(SessionName1,SessionName2);
    fname = strcat(sData1.sessionInfo.fileID,'-roi',num2str(i),'_PosTuning');
    savefig(fullfile(SavePath,fname));
    saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));
 end
close all

% plot mean
PlaceCellComp.MeanPosTuningBothS1 = mean(MaoPC1.PosTuning(PlaceCellComp.BothPlaceCells.ROIs,:),1);
PlaceCellComp.MeanPosTuningBothS2 = mean(MaoPC2.PosTuning(PlaceCellComp.BothPlaceCells.ROIs,:),1);
figure();
Ymax = max(max(PlaceCellComp.MeanPosTuningBothS1),max(PlaceCellComp.MeanPosTuningBothS2))+0.0000001;
plot(Xaxis,PlaceCellComp.MeanPosTuningBothS1,'LineWidth',2); hold on; % pos tuning curve
plot(Xaxis,PlaceCellComp.MeanPosTuningBothS2,'LineWidth',2); hold on; % pos tuning curve
axis([0 160 0 Ymax]); % ceil(Ymax)
title(strcat(sData1.sessionInfo.fileID,' Pos. tuning MEAN both'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
legend(SessionName1,SessionName2);
fname = strcat(sData1.sessionInfo.fileID,'-roi',num2str(i),'_PosTuning-MEAN(both)');
savefig(fullfile(SavePath,fname));
saveas(gcf,(fullfile(SavePath,[fname '.jpg'])));

%save(fullfile(sData1.sessionInfo.savePath,strcat(sData1.sessionInfo.fileID,'_PlaceCellComp.mat')),'PlaceCellComp'); C:\MATLAB\SAVE\
save(fullfile('C:\MATLAB\SAVE\',strcat(sData1.sessionInfo.fileID,'_PlaceCellComp.mat')),'PlaceCellComp');

%{
%%% cue vs non-cue region analysis
CueCellsPre1 = PlaceCellComp12.BothPlaceCells.PosTuningPeak;
CueCellsPre2 = CueCellsPre1(CueCellsPre1(:,2)>8,1:4);
CueCellsBoth = CueCellsPre2(CueCellsPre2(:,2)<45,1:4);
NonCueCellsPre1 = CueCellsPre1(CueCellsPre1(:,2)<8,1:4);
NonCueCellsPre2 = CueCellsPre1(CueCellsPre1(:,2)>45,1:4); 
NonCueCellsBoth = vertcat(NonCueCellsPre1,NonCueCellsPre2);

ShiftInCueCells = mean(CueCellsBoth(:,4));
ShiftInNonCueCells = mean(NonCueCellsBoth(:,4));
%}

end