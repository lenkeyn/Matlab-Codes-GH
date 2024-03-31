%%% NORI ANALYSIS PIPLINE - SDATA ERA

%%%%% ANALYSIS OF IMAGED CA-TRANSIENTS AND BEHAVIOR
                            
%%%%% MOTION CORRECTION OF CA-SESSIONS
% STEP #1 : SET THE FOLDER IN WHICH YOUR SCISCAN RAW DATA CAN BE FOUND (IT WILL CHOOSE AUTOMATICALLY THE CORRECT FILE) 
[~,filePath,~] = uigetfile('*.raw');
filePath = fileparts(filePath);

% STEP #2 : LOAD THE IMAGES AND METADATA INTO MATLAB
meta2P = getSciScanMetaData(filePath);

% STEP #3 : RIGID MOTION CORRECTION PERFORMED AUTOMATICALLY, BUT SHIFTS SHOULD BE SAVED TO THE DATA FOLDER TOGETHER WITH ALIGNED IMAGES
% parameters to be set in registerImages.m, examples in demo.m.
% options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',50,'init_batch',200);
% options_nonrigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'grid_size',[32,32],'mot_uf',4,'bin_width',200,'max_shift',15,'max_dev',3,'us_fac',50,'init_batch',200);
% Eivind suggestion for non-rigid: bin_width 10-20, init_batch 50-100,
% grid_ size usually 128, so 128/64/32 is OK (but smaller is slower).
 
registerImagesNori(filePath); 
registerImagesNoNRigid(filePath);
registerImagesNRNori(filePath);

% STEP #4 : CREATE ROIS IN roimanager_lite, SAVE ROIs, SAVE SIGNAL
roimanager; 
  

%%%%% CONVERSION AND ANALYSIS OF LABVIEW DATA
%%% CREATE mouseInfo struct into C\MATLAB\MOUSEINFO folder (once for each mouse)
% Use: CreateMouseSheet_empty.mat (not a function)
% sData = CalcBehavOldMicheleConvert(nBins); % Michele's LVdata use


%%% Create Session Info, DAQData Info, Behavioral Info:
IsMoreStim = 0; % were there different opto-stim protocols?
RedLaserMax = 11; % 100% laser stim mW/mm2. Big beam: 10-11 mm2, and full intensity usually set to 10 mW, small beam 1mm2,  
IsOpto = 1;
[sData,filePath,fileName] = CreateSessionInfo_LV91(IsOpto); %Log text file recorded during experiment will be opened and automatically looks for mouseinfo
%%% CONVERT LABVIEW RAW DATA INTO TDMS WHICH CAN BE USED BY MATLAB (Lick, PD, Water, Frame, Counter)
sData = loadTDMSdataNori(sData,filePath,fileName);
%%% calculate behavioral variables
%sData.daqdata.frameSignal(:,1) = 0;

nBins = 80; % number of bins
%IsOptoSorted = 1; % was it an optical stimulated session? 0:no, 1:yes
OptoStimLimitMs = 15000; % in ms
DiscardCmBeginningCm = 10;
OptoSensitivity = 5; % set 5 for sinus fT (also if stim goes to the next lap, until reward), set to 100 to short stimulations - then very sensitive
sData = CalcBehav2(sData,nBins,IsOpto,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm);   % sData.behavior.opto.IsOptoSession = 1;

%%% Comparing First - Second half - only behavior
%IsOptoSorted = 1; % if all trials was laser on or off
LapsTested = 20;
sData = BehavFirstSecondHalf(sData,LapsTested); 
close all

%optoMoreProts = 0;
%DiscardBeginningCm = 10;
%sData = OptoTrialSortingMoreStim2(sData,DiscardBeginningCm); % if more than one optical stimulation was used during session. Failed trials labelled with -1, dark grey on plot.

%%% ANALYSING INDIVIDUAL CALCIUM IMAGING SESSION:
% after extracting data from ROI manager you will have to open:
% _signals.mat_cia_deconvolved , _signals.mat_cia_denoised , _signals.mat_dff , _signals.mat_cia_spikes
% Use data extracted from RoiManager: open dff data, sData
VelMin = 0.1; %SET!!!! interneurons = 0 , PC= 0.1
pmtGainGreen = 25;%!!!
waveLength = 920;
laserPower = 85; % mW
fovCoord = [774 671];          % AP, ML coordinated from center of window (0 -500 means imaging at -2.2 mm left hemisphere, ML +right, - left hem) left right changed 2019.04.03.
IsDeconv = 1; % 2: do the deconvolution, 1: use ROI manager deconvolution, 0 : do not use deconvolution
gol = 5; % gol#0 original reward, gol#2 50 cm forward, between hot glues, gol#3 50 cm before original after Velcro, gol = 10 no cues
sData = calcCaDataNori2(sData,VelMin,laserPower,waveLength,fovCoord,pmtGainGreen,IsDeconv,gol); %calcCaDataNori %calcCaDataNoriShortVIP
clearvars IsOpto BinNu VelMin  laserPower waveLength fovCoord
close all
 
%Sorting Optodata and plotting heatmap and mean act for opto-on/off/after trials:
%sData =  CaDataOptoSortingdFFsData2(sData); 
FigVisible = 'off';
sData = CaDataOptoSortingdFFsData2(sData,FigVisible); 
%sData  = CaDataOptoSortingDeconvData2(sData); % no smoothing
%sData  = CaDataOptoSortingCumSpikesData2(sData); % no smoothing
%{
% if more than one optical stimulation was used during session
FigVisible = 'off'; datatype = 0;
ProtocolNames = string({'none','stim 1 mW/mm2','stim 0.3 mW/mm2'}); % check protocol in sData.stimProtocols.protocol
%'stim 10 mW','stim 6mW','stim 3 mW ','stim 1 mW','stim 14-84 cm','stim 86-156 cm'
sData  = CaDataOptoSortingdFFsDataMoreProts2(sData,datatype,FigVisible,ProtocolNames);  

%%% Place cell anal for non-opto sessions
datatype = 0; % use dff data , datatype = 1 for deconvolved data
sData = placeCellMaosData2(sData,datatype); % final params: activityTreshold = 0.5; MinPlaceFieldSize = 2; MaxPlaceFieldSize = 100; InOutRatio = 2; ReliabiliyIndex = 0.3;
sData = LandmarkCellDetection(sData); % detects cells which respond to both cues (e.g. both hot glue spikes or both velcro), detects place cells with only one place field as well
%}

%%% Place cells in optically stimulated sessions:
%InOutRatio = 2; %ReliabiliyIndex = 0.3;
datatype = 0; % use dff data , datatype = 1 for deconvolved data
FigVisible = 'off';
sData = placecellMaosDataOptoOn2(sData,datatype,FigVisible); %2024.01.01. input: savePath %sData = placecellMaosDataOptoOn2(sData,datatype); placecellMaosDataOptoOn3 is a previous version
sData = LandmarkCellDetectionOptoOnDFFVIPpaper3(sData); %2024.01.01.  detects cells which respond to both cues (e.g. both hot glue spikes or both velcro), detects place cells with only one place field as well

%{
nProtOpto = 2;
for i = 1:1:nProtOpto+1
    Protocol = i; %Protocol 1: ctr, protocol 2-3-4 opto, Protocol 5: after opto
    datatype = 0; InOutRatio = 2; ReliabiliyIndex = 0.3; FigVisible = 'off';
    sData = placecellMaosDataOptoOnMoreProt2(sData,datatype,InOutRatio,ReliabiliyIndex,Protocol,FigVisible);
end
Protocol = 5; % Protocol 5: after opto
datatype = 0; InOutRatio = 2; ReliabiliyIndex = 0.3; FigVisible = 'off';
sData = placecellMaosDataOptoOnMoreProt2(sData,datatype,InOutRatio,ReliabiliyIndex,Protocol,FigVisible);


% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

%%% compare the same place cells in different sessions
plotmeanROIActAllSessions2(sData,IsOptoSession,type) %%% have to set sctipt!
%}

%%% compare the same place cells in different types of protocol within one session (opto-on, opto-off)
type = 0; % 0 dff % type = 1; % deconv
GaussFilter = 5;
sData = placecellComparison4_1stim(sData,type,GaussFilter); %2024.01.01.

%type = 0;
%sData = placecellComparison3_Morestim(sData,type);

%%% opto effect 
%sData = effectOnPCCtrVIPpaper2(sData); 
sData = effectOnPCCtrOnePFVIPpaper3(sData); %2024.01.01.
sData = effectOnPCLandmarkCellsVIPpaper3(sData); %2024.01.01.
sData = effectOnPCBeforeRewardVIPpaper3(sData); %2024.01.01.
sData = EffectPCinCtr_VIPPAPER(sData); %2024.01.01. effect on cell level 

%%% new place field analysis
%newPlaceFields(sData)
sData = newPlaceFields4_PFStartBin(sData);
% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

%%% additional new codes:
gainModulationInOutRatioPCinCtr(sData);
%gainModulationInOutRatioLandmarkAndOnePFPC(sData);
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

%%% for VIP data, recording the same cells over days - NOT UPDATED
%ROIsDiscard = sData.imdata.ROIsToDiscard; % ROIs = []; analyse all ROIs
sData = placeCellMaosDataROIs(sData,type,InOutRatio,ReliabiliyIndex,sData.imdata.ROIsToDiscard);
close all;
