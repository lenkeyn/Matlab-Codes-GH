function [BinLickMatrix,MeanLickBin,LickCmMatrix,MeanLickCm] =  plotHeatBinLicksPlusBin(LVDATA,AllBinNu) 

%PARAMETERS TO BE SET:
% BinSize = LVDATA.BinSize; % used for matrix
% SampleNu = LVDATA.FRNu; % used for matrix 
TRNu = LVDATA.TRNu;
BinNu = AllBinNu;
EnterIntoBin = LVDATA.EnterIntoBinSampleIndPlusBins;
SampleSpentInBin = LVDATA.SampleSpentInBinPlusBins;
% BinNu = floor(50*pi/BinSize);  % Circumference = 50*pi; %wheel circumference in cm, diameter is 50 cm

% READ DATA FROM FILE
% Distance = LVDATA.AbsDistDSMonIncr; % used for matrix
% Licking = TDMSDATA.Licksignal(LVDATA.FrameStart); % Eivind downsapled in this way
Temp = 0; % insert an extra zero to the beginnning
Licking = vertcat(Temp,LVDATA.LickDS); % check if there is lick within at least half the duration if the frame scan, if yes, consider lick (=1)

% set trial number for plotting, search last full data trial
for i = 1:1:TRNu
    if any(isnan(EnterIntoBin(i,:))) || any(EnterIntoBin(i,:)==0)
        TRNuPlot = i-1; 
        break
    else
        TRNuPlot = i;
    end
end


BinLickMatrix = NaN(TRNuPlot,BinNu);
% calculate the number of licks during each bin (lick/cm)
for i = 1:1:TRNuPlot  % rows are trials
    for j = 1:1:BinNu  % columns (distance bins)  
        LR = 0; % set lickrate to zero before counting in each bin
        for m = 1:1:SampleSpentInBin(i,j) % Note: if the recorging starts with lick, it will not detect that first lick
            if Licking(EnterIntoBin(i,j) + m - 2) == 0 && Licking(EnterIntoBin(i,j) + m-1) == 1   % search the first timepoint (sample) when licking starts, zero changes to one, this will be a lick event
                LR = LR + 1;  % if there is a lick event, increase lick rate with one
            end
        end
        BinLickMatrix(i,j) = LR; % put the cumulative lick number in a bin into BinLick matrix
    end
end
MeanLickBin = nanmean(BinLickMatrix,1);

LickCmMatrix = BinLickMatrix./LVDATA.BinSize;
MeanLickCm = nanmean(LickCmMatrix,1);
%Cmax = ceil(max(MeanLickCm)*2);

%PLOT FIGURE
figure('Color','white'); 
imagesc(1:LVDATA.BinSize:(AllBinNu*LVDATA.BinSize),1:TRNuPlot,(LickCmMatrix)) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
j = colorbar;
colormap(jet);
j.Label.String = 'Licking (licks/cm)';
j.Label.FontSize = 11;
j.TickDirection = 'out'; 
caxis([0 2]); %set limits for color plot, below 1st black, above 2nd white
hold on;
line([157 157],[0 TRNuPlot],'Color','white');
hold on
line([163 163],[0 TRNuPlot],'Color','white');

xlabel('Position on wheel (cm)');
ax = gca;
ax.TickDir = 'out';
%xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
%xticks([0,25,50,75,100,125,150]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
if TRNuPlot < 50
    yticklabels = 0:5:TRNuPlot;
elseif TRNuPlot >= 50 && TRNuPlot < 200
    yticklabels = 0:10:TRNuPlot;
else
    yticklabels = 0:20:TRNuPlot;
end
yticks = linspace(1, TRNuPlot, numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(LVDATA.FileID);

end