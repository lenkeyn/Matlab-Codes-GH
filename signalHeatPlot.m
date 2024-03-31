function signalBinned = signalHeatPlot(signal,sData,VelMin,FigVisible)

% create binned data

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:sData.imdata.nROIs
    signalBinned{i} = NaN(sData.behavior.wheelLapImaging-1,sData.behavior.meta.nBins); 
end
% load matrices with data in trial and bins
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:));
% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBin Matrix
SampleInBinVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end
% CREATE MATRIX for Velo-limited data based on velocity limit (discarded data will be represented by -1)
SampleInBinLim = sData.behavior.binning.samplesInBinIndex; % CADATA.SampleInBinLim = CADATA.SampleInBin;
% substitute indexes with low velo data with -1 : 
SampleInBinLim(SampleInBinVelo < VelMin) = -1;
for i = 1:1:sData.imdata.nROIs   
    for j = 1:1:sData.behavior.wheelLapImaging-1
        for k = 1:1:sData.behavior.meta.nBins
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBin = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            dFFInBin = NaN(MaxSamplesinBin,1); % temporary array to calcuate mean Ca-value in a bin during time spent in a bin (> velo lim)
            for m = 1:1:SpentInBin
                if SampleInd+m > sData.imdata.nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  dFFInBin(m) = signal(i,(SampleInd+m-1)); % dFF
                end
            end
            signalBinned{i}(j,k) = nanmean(dFFInBin); % mean of Ca data within a bin witihn a trial
        end
    end
end

%{
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-signal');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-signal'); 
MeanRoiActSignal = NaN(sData.imdata.nROIs,sData.behavior.meta.nBins);
%FigVisible = 'on';
for roi = 1:1:sData.imdata.nROIs
    MeanRoiActSignal(roi,1:sData.behavior.meta.nBins)= nanmean(signalBinned{roi},1);
    if(any(isnan(MeanRoiActSignal(roi,1:sData.behavior.meta.nBins))))
        continue 
    end
    plotHeatBinCa(signalBinned{roi},sData.sessionInfo.fileID,roi,'activity',sData.behavior.meta.binSize,sData.behavior.meta.nBins,sData.behavior.wheelLapImaging-1,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 sData.behavior.wheelLapImaging-1],'Color','white','LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 sData.behavior.wheelLapImaging-1],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 sData.behavior.wheelLapImaging-1],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 sData.behavior.wheelLapImaging-1],'Color','white','LineStyle','--','LineWidth',2); 
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-signal');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

    figure('Color','white','visible',FigVisible')
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*sData.behavior.meta.nBins;
    Ymax = (max(MeanRoiActSignal(roi,:)))+0.0001;
    plot(Xaxis,MeanRoiActSignal(roi,1:sData.behavior.meta.nBins),'LineWidth',2)
    line([sData.imdata.cues.C1A sData.imdata.cues.C1A],[0 sData.behavior.wheelLapImaging-1],'Color','black','LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C2A sData.imdata.cues.C2A],[0 sData.behavior.wheelLapImaging-1],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C3A sData.imdata.cues.C3A],[0 sData.behavior.wheelLapImaging-1],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); 
    line([sData.imdata.cues.C4A sData.imdata.cues.C4A],[0 sData.behavior.wheelLapImaging-1],'Color','black','LineStyle','--','LineWidth',2);   axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,25,50,75,100,125,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning-signal');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
end
%}
%sData.imdata.binned.MeanRoiAct_signal = MeanRoiActSignal;
close all;

    
end