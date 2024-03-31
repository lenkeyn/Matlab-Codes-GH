function sData = CaDataonlySR(sData)

%%% TO BE SET:  
% VelMin = 0.01; % minimum velocity. Below this, Ca-activity will be discarded
% IsDeconv = 1; % 2: do the deconvolution, 1: use ROI manager deconvolution, 0 : do not use deconvolution
% gol = 5; % gol#0 original reward, gol#2 50 cm forward, between hot glues, gol#3 50 cm before original after Velcro, gol = 10 no cues


%%% SET PARAMETERS 
VelMin = 0;
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging;
%FrameRate = sData.behavior.meta.imagingSamplingRate;
[nROIs, nSamplesPre] = size(sData.imdata.roiSignals(2).dff);
nSamples = nSamplesPre -5;

% saving folder
mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'RoiActBinned-SR');

% calculate slow removed transients, filtering Window in sec, percentile is for calculating baseline (ususlly 10-20%)
Window = 10; % based on m8056 VIP-Gcamp recording I find 10 sec good
Percentile = 20; % 20%
sData = SlowRemoved(sData,Window,Percentile);


%%% PLOT TRANSIENTS
%savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
% calculate lowpassed data for visulization
% LIGHT filtering to see better fast transients
dff_lowpassLightSR = NaN(nROIs,nSamples);
d = fdesign.lowpass('Fp,Fst,Ap,Ast',0.1,1,1,60); % Potential changes for stronger filtering: decrease Fp (0.01) and Fst (0.08). original stting: 0.1,1,1,60
Hd = design(d,'cheby1');   % or use 'equiripple' filter to have less amplitude filtering and more noise
for i = 1:1:nROIs
    %dff_lowpassLight(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff(i,1:nSamples)); 
    dff_lowpassLightSR(i,:) = filter(Hd,sData.imdata.roiSignals(2).dff_slowRemoved(i,1:nSamples));
end
%sData.imdata.roiSignals(2).dff_LPlight = single(dff_lowpassLight);
sData.imdata.roiSignals(2).dff_slowRemoved_LPlight = single(dff_lowpassLightSR);

% plot slightly filtered dff_slowremoved transients
%nROIOnFig = round((nROIs/3),-1)-1; % how many ROIs to see on the fig (best to set to X if rem(X/5)=1 , e.g. 19, 49, ...)
% plot 10 ROIs to one fig
SavePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiAct');
nROIOnFig = 10; 
plotdata = dff_lowpassLightSR(:,1:nSamples);
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIOnFig,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,SavePath,'SR-NormROIsdFF-LPlight');
%plot all rois into one fig
plotMultipleROIdFFNormAbsDist(plotdata,sData.behavior.meta.imagingSamplingRate,nROIs,sData.behavior.wheelPosDsMonIncr,sData.behavior.lickDs,SavePath,'SR-ALL-NormROIsdFF-LPlight');  %sData.behavior.opto.LightOnSignalDS
close all;

% save temp
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');


%%% CREATE A NEW MATRIX CONTAINING REAL (all) VELOCITY DATA , same size as SampleInBin Matrix
SampleInBinVelo = NaN(size(sData.behavior.binning.samplesInBinIndex));
% load matrix with velo data
for i = 1:1:(size(sData.behavior.binning.samplesInBinIndex,1))-1   % -1    find in which bin this sample is 
    SampleInBinVelo(i,sData.behavior.binning.samplesInBinIndex(i,:)==i) = sData.behavior.runSpeedDs(i); % replace sample index to velocity 
end

% CREATE MATRIX for Velo-limited data based on velocity limit (discarded data will be represented by -1)
SampleInBinLim = sData.behavior.binning.samplesInBinIndex; % CADATA.SampleInBinLim = CADATA.SampleInBin;
% substitute indexes with low velo data with -1 : 
SampleInBinLim(SampleInBinVelo < VelMin) = -1; % CADATA.SampleInBinLim(CADATA.SampleInBinVelo < VelMin) = -1;

% PUT SAMPLE INDEXES (Velo limited) into multi-dimension matrix. 1st dim: trials, then columns are bins, rows contain samples in that trial and bin. 
MaxSamplesinBin = max(sData.behavior.binning.samplesSpentInBin(:)); % how many rows needed maximum in matrix (what was max time (samples) spend in a bin)   , MaxSamplesinBin = max(CADATA.SampleSpentInBin(:));
for i = 1:1:nTrials
    MMSamplesLim{i} = NaN(MaxSamplesinBin,nBin); % for samples, creates i times a row-col matrices (row: max number of samples spent in a bin, column: Bin number)
    MMVelo{i} = NaN(MaxSamplesinBin,nBin); % same for all Velo data (wo limitation)
    MMVeloLim{i} = NaN(MaxSamplesinBin,nBin); % same for limited Velo data
end
% FILL UP MATRICES + Calculate average velocity within a bin and put into  CADATA.VeloInBin and CADATA.VeloLimInBin matrix
sData.imdata.VeloInBin = NaN(nTrials-1,nBin); % mean real velo
sData.imdata.VeloLimInBin = NaN(nTrials-1,nBin); % mean limited velo
for i = 1:1:nTrials-1
   for j = 1:1:nBin
      SampleInd = sData.behavior.binning.enterIntoBinIndex(i,j); % get sample index when enter into a bin  % SampleInd = CADATA.EnterIntoBinSampleInd(i,j); 
      if isnan(SampleInd)
         break 
      end
      SpentInBin = sData.behavior.binning.leaveBinIndex(i,j) - sData.behavior.binning.enterIntoBinIndex(i,j)+ 1 ; % get how many samples spent in that bin  % SpentInBin = CADATA.SampleSpentInBin(i,j);
      for k = 1:1:SpentInBin
          if SampleInBinLim(SampleInd+k-1,j) > 0 % limited values were set to -1, I do not want to contain them
             MMSamplesLim{i}(k,j) = SampleInd + k - 1; % limited samples
             MMVeloLim{i}(k,j) = SampleInBinVelo((SampleInd+k-1),j); % use only limited velo data  % CADATA.MMVeloLim{i}(k,j) = CADATA.SampleInBinVelo((SampleInd+k-1),j);
          end
          MMVelo{i}(k,j) = SampleInBinVelo((SampleInd+k-1),j); % use all velo data  % CADATA.MMVelo{i}(k,j) = CADATA.SampleInBinVelo((SampleInd+k-1),j); 
      end
      sData.imdata.VeloLimInBin(i,j) = nanmean(MMVeloLim{i}(:,j)); % mean of limited velo
      sData.imdata.VeloInBin(i,j) = nanmean(MMVelo{i}(:,j));  % mean of real velo
   end
end

% CREATE A NEW MATRIX FOR Ca-data WHERE ALL DATA IN A BIN WITHIN A TRIAL IS REPRESENTED BY ONE CELL/NUMBER. CA-DFF DATA WILL BE AVERAGED WITHIN ONE BIN/ONE TRIAL (USING DATAPOINTS WHERE SPEED IS HIGHER THAN SET)
for i = 1:1:nROIs
    sData.imdata.binned.RoidFFSR{i} = NaN(nTrials,nBin);
end
 % load matrices with data in trial and bins
for i = 1:1:nROIs   
    for j = 1:1:nTrials-1
        for k = 1:1:nBin
            SampleInd = sData.behavior.binning.enterIntoBinIndex(j,k); % get sample index when enter into a bin
            if isnan(SampleInd) % if recording ends in LV stop calculation
               break 
            end
            SpentInBin = sData.behavior.binning.leaveBinIndex(j,k) - sData.behavior.binning.enterIntoBinIndex(j,k) + 1; % get how many samples spent in that bin
            dFFSRInBin = NaN(MaxSamplesinBin,1);
            for m = 1:1:SpentInBin
                if SampleInd+m > nSamples % if recording ends in SciScan stop calculation (sometimes not the same size data in SciScan and LV)
                    break
                end
                if SampleInBinLim(SampleInd+m-1,k) > 0 % limited values were set to -1, I do not want to contain them
                  dFFSRInBin(m) = sData.imdata.roiSignals(2).dff_slowRemoved(i,(SampleInd+m-1));
                end
            end
            sData.imdata.binned.RoidFFSR{i}(j,k) = nanmean(dFFSRInBin); % mean of Ca data within a bin witihn a trial
        end
    end
end
%}

nTrials = sData.behavior.wheelLapImaging;
nBins = sData.behavior.meta.nBins;
nROIs = sData.imdata.nROIs;
% plotdata: sData.imdata.binned.RoidFF_SR_LP; sData.imdata.binned.RoidFF_SR; sData.imdata.binned.RoidFF

C1A = 26; C1B = 26;
C2A = 46; C2B = 46;
C3A = 66; C3B = 66;
C4A = 86; C4B = 86;
roiStart = 1;
roiEnd = nROIs;
%%%Plot dFF_slowremoved signals
FigVisible = 'off';
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\RoiActBinned-SR'); % '\Imaging\RoiActBinned_SRLP'
MeanRoiAct = NaN(nROIs,nBins);
for roi = roiStart:1:roiEnd %roiStart
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFFSR{roi},1);

    if(any(isnan(MeanRoiAct(roi,1:nBins))))
       continue 
    end
    % sData.imdata.binned.RoidFF_SR
    plotHeatBinCa(sData.imdata.binned.RoidFFSR{roi},sData.sessionInfo.fileID,roi,'dF/F(SR)',sData.behavior.meta.binSize,sData.behavior.meta.nBins,nTrials,FigVisible); %plotHeatBinCa(data,fileID,roi,ylab,BinSize,nTrials,figure visibility)
    caxis([0 inf]); hold on;
    line([C1A C1A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 nTrials],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 nTrials],'Color','white','LineStyle','--','LineWidth',2); hold on;
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-dff');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
  
    figure('Color','white','visible',FigVisible'); % 'visible','off'
    Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
    % sData.imdata.binned.RoidFF_SR_LP
    MeanRoiAct(roi,1:nBins)= nanmean(sData.imdata.binned.RoidFFSR{roi},1);
    Ymax = (max(MeanRoiAct(roi,:)))*1.1;
    plot(Xaxis,MeanRoiAct(roi,1:nBins),'LineWidth',2)
    line([C1A C1A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    line([C2A C2A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C3A C3A],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[0 Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
    line([C4A C4A],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[0 Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
    axis([0 160 0 Ymax]); % ceil(Ymax)
    title(strcat(sData.sessionInfo.fileID,' ROI #',num2str(roi)));
    xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
    ylabel('Position tuning of activity');
    fname = strcat(sData.sessionInfo.fileID,'-roi',num2str(roi),'-pos-tuning');
    savefig(fullfile(savePath,fname));
    saveas(gcf,(fullfile(savePath,[fname '.jpg'])));
    %}
end
close all;
sData.imdata.binned.MeanRoiActSR = MeanRoiAct;


%%% PLot Mean position tuning of all ROIs slow removed transeints
figure('Color','white');
Xaxis = sData.behavior.meta.binSize:sData.behavior.meta.binSize:sData.behavior.meta.binSize*nBins;
sData.imdata.binned.MeanMeanRoiActSR = nanmean(sData.imdata.binned.MeanRoiActSR);
Ymax = (max(sData.imdata.binned.MeanMeanRoiActSR(1,:)))*1.1;
Ymin = (min(sData.imdata.binned.MeanMeanRoiActSR(1,:)))*0.9;
plot(Xaxis,sData.imdata.binned.MeanMeanRoiActSR(1,1:nBins),'LineWidth',2)
line([C1A C1A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C1B C1B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
line([C2A C2A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C2B C2B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C3A C3A],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on; line([C3B C3B],[Ymin Ymax],'Color',[1 0.3 0.5],'LineStyle','--','LineWidth',2); hold on;
line([C4A C4A],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on; line([C4B C4B],[Ymin Ymax],'Color','black','LineStyle','--','LineWidth',2); hold on;
axis([0 160 Ymin Ymax]); % ceil(Ymax)
title(strcat(sData.sessionInfo.fileID,' Mean of all ROIs SR'));
xlabel('Position on Wheel (cm)'); ax = gca; ax.TickDir = 'out'; xticks([0,50,100,150]);
ylabel('Position tuning of activity');
fname = strcat(sData.sessionInfo.fileID,'AllRois-pos-tuning-SR');
savefig(fullfile(savePath,fname));
saveas(gcf,(fullfile(savePath,[fname '.jpg'])));

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end

