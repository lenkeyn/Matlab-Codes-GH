function [VeloMatrix,MeanVeloBin] = plotHeatBinVelo(LVDATA,Ymax) 

%SET FOR EACH PLOT!
%BinSize = 3; % stepsize in cm for binning

%PARAMETERS TO BE SET:
%BinSize = LVDATA.BinSize; % used for matrix
%SampleNu = LVDATA.FRNu; % used for matrix 
TRNu = LVDATA.TRNu;
BinNu = LVDATA.BinNu;
EnterIntoBin = LVDATA.EnterIntoBinSampleIndPlusBins;
SampleSpentInBin = LVDATA.SampleSpentInBinPlusBins;
%BinNu = floor(50*pi/BinSize);  % Circumference = 50*pi; %wheel circumference in cm, diameter is 50 cm

% READ DATA FROM FILE
%Distance = LVDATA.AbsDistDSMonIncr; % used for matrix
Velocity = LVDATA.VelDS;

% calculating sample in bin matrix. Now it is calculating during LVDATA calculation
%{
SampleInBin = NaN(SampleNu,BinNu);    
for j = 1:1:(BinNu-1) % first search all data for the specified distance bin (indexes for entering a new bin can be calculated based on AbsDistDS data, j is distance cm, bin on the wheel (e.g. 0-1cm)
    for i = 1:1:SampleNu % scan through the whole dataset and sort into bin j, between trials, there are NaNs on the column
       if  Distance(i) >= (j-1)*BinSize && Distance(i) < j*BinSize  %check ALL distanceDS datapoint if it belongs to the specified bin, if yes, continue
           SampleInBin(i,j) = i; % put sample index into the bin column. 
       end
    end 
end
if j == BinNu-1 % last Bin, which is not a full bin or bigger, collect together what is before lapstart into last bin
   for i = 1:1:SampleNu
      if  Distance(i) > j*BinSize 
          SampleInBin(i,j+1) = i;          
      end
   end  
end


% Make 2 matrices. CADATA.EnterIntoBinSampleInd: when the animal enter into a bin (sample ind, row: trial, col: bin); SampleSpentInBin: how many samples spend the aimal in this bin. At the end of recording there are NaNs in the matrices 
SampleInBinIsNaN1 = isnan(SampleInBin);
SampleInBinIsNaN2 = diff(SampleInBinIsNaN1,1)==-1; % matrix shows first sample in a bin. It is 1 and other is 0. (row: trial, column: bin)
SampleInBinIsNaN3 = diff(SampleInBinIsNaN1,1)==1; % matrix when last sample spent in a given bin, set to 1 and other is 0.
LastSampleFix = find(SampleInBinIsNaN1(SampleNu,:) == 0);  % I have to drop the last sample in order the code to function prperly
SampleInBinIsNaN3(SampleNu-1,LastSampleFix) = 1 ; % drop last sample 
EnterIntoBinSampleInd = NaN(TRNu,BinNu); % sample index when animal enter into given bin given trial
LeaveBinSampleInd = NaN(TRNu,BinNu); % how many samples spent in a given bin in a given trial
for i = 1:1:BinNu % I need to do it in a complicated way using temporary arrays because if the last trial did not go to the end it gave error (mismatch in TRNu and bin-start at later bins)
    TempArray1 = zeros(TRNu,1);
    TempArray1 = find(SampleInBinIsNaN2(:,i)==1)+1; % first sample spend in a given bin, given trial
    if numel(TempArray1)<TRNu
        TempArray1(TRNu) = NaN;
    end
    EnterIntoBinSampleInd(:,i) = TempArray1;
    TempArray1 = zeros(TRNu,1);
    TempArray1 = find(SampleInBinIsNaN3(:,i)==1); % last sample spend in a given bin, given trial
    if numel(TempArray1)<TRNu
        TempArray1(TRNu) = NaN;
    end
    LeaveBinSampleInd(:,i) = TempArray1; % last sample spend in a given bin, given trial
end
SampleSpentInBin = LeaveBinSampleInd - EnterIntoBinSampleInd + 1; % spent in bin = last sample - first sample +1

clearvars LastSampleFix SampleInBinIsNaN1 SampleInBinIsNaN2 SampleInBinIsNaN3 TempArray1 row col; % cleaning up workspace
%}

% set trial number for plotting, search last full data trial
for i = 1:1:TRNu
    if isnan(EnterIntoBin(i,size(EnterIntoBin,2))) || EnterIntoBin(i,size(EnterIntoBin,2))==0
        TRNuPlot = i-1; 
        break
    else
        TRNuPlot = TRNu;
    end
end
%{
VeloMatrix = NaN(TRNuPlot,BinNu);
% calculate the velo during each bin 
for i = 1:1:TRNuPlot  % rows are trials
    for j = 1:1:BinNu  % columns (distance bins)  
        VeloMatrix(i,j) = mean(Velocity(EnterIntoBin(i,j):LeaveBin(i,j)));
    end
end
MeanVeloBin = nanmean(VeloMatrix,1);
%}


VeloMatrix = NaN(TRNuPlot,size(EnterIntoBin,2));
for i = 1:1:TRNuPlot  % rows are trials
    for j = 1:1:size(EnterIntoBin,2)  % columns (distance bins)  
      VeloMatrix(i,j) = mean(Velocity(EnterIntoBin(i,j):(EnterIntoBin(i,j)+SampleSpentInBin(i,j)-1)));
    end
end
MeanVeloBin = nanmean(VeloMatrix,1);


%PLOT FIGURE
figure('Color','white'); 
imagesc(1:160,1:(TRNuPlot),(VeloMatrix(1:TRNuPlot,1:BinNu))) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
j = colorbar;
colormap(jet);
j.Label.String = 'Velocity (cm/s)';
j.Label.FontSize = 11;
j.TickDirection = 'out'; 
caxis([0 Ymax]); %set limits for color plot, below 1st black, above 2nd white

xlabel('Position on wheel (cm)');
ax = gca;
ax.TickDir = 'out';
%xticklabels = 0:20:160; %cannot set the axis to show 157 as the end%
xticks([0,25,50,75,100,125,150]);
set(gca, 'XTick', xticks, 'XTickLabel', xticklabels)

ylabel('Trials');
if TRNuPlot < 50
    yticklabels = 0:5:TRNuPlot;
elseif TRNuPlot >= 50 && TRNuPlot < 200
    yticklabels = 0:10:TRNuPlot;
else
    yticklabels = 0:20:TRNuPlot;
end
yticks = linspace(1, (TRNuPlot), numel(yticklabels));
set(gca, 'YTick', yticks, 'YTickLabel', yticklabels)
title(LVDATA.FileID);


end