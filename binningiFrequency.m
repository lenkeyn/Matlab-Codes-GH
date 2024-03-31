function sData = binningiFrequency(sData) % signal: which signal to be binned: deconv, iFreq

% calculate calcium transient instantenous frequency
%sData = ifreqDeconv(sData);

% calculate binned data of 
%signal = sData.imdata.roiSignals(2).iFrequency; 
str1 = 'iFreq';

%mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),strcat('RoiActBinned-',str1));
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-',str1); 

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
nSamples = sData.imdata.meta.nFrames;
VelMin = 0.1;
 

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdata.binned.RoiiFrequency{i} = NaN(nTrials,nBins); 
end
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:));
%%% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBins Matrix
SampleInBinVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end
SampleInBinsLim = sData.behavior.binning.samplesInBinIndex; 
SampleInBinsLim(SampleInBinVelo < VelMin) = -1; 

 % load matrices with data in trial and bins
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:nBins
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBins = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            iFreqInBins = NaN(MaxSamplesinBin,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
            for m = 1:1:SpentInBins
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinsLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  iFreqInBins(m) = sData.imdata.roiSignals(2).iFrequency(i,(SampleInd+m-1)); % dFF
                end
            end
            sData.imdata.binned.RoiiFrequency{i}(j,k) = nanmean(iFreqInBins); % mean of Ca data within a bin witihn a trial
        end
    end
end

 
roiStart = 1;
roiEnd = nROIs;
% Bin data
FigVisible = 'off';
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoiiFrequency{roi},1);
    %{
    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    plotHeatBinCa(sData.imdata.binned.RoiiFreq{roi},sData.sessionInfo.fileID,roi,'instantanenous frequency of activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-iFreq');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    % sData.imdata.binned.RoidFF_SR_LP
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoiiFrequency{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1+0.0001;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-iFreq');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
close all;
sData.imdata.binned.MeanRoiAct_iFreq = MeanRoiAct;

end