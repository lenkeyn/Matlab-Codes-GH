%%% NORI ANALYSIS PIPLINE - SDATA ERA

%%%%% ANALYSIS OF IMAGED CA-TRANSIENTS AND BEHAVIOR

%%%%% MOTION CORRECTION OF CA-SESSIONS
% STEP #1 : SET THE FOLDER IN WHICH YOUR SCISCAN RAW DATA CAN BE FOUND (IT WILL CHOOSE AUTOMATICALLY THE CORRECT FILE) 
[~,filePath,~] = uigetfile('*.raw');

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
sData = CalcBehavOldMicheleConvert(BinNu); % Michele's LVdata use

%%% Create Session Info, DAQData Info, Behavioral Info:
IsMoreStim = 1; % were there different opto-stim protocols?
RedLaserMax = 10; % 100% laser stim mW/mm2. Big beam: 10-11 mm2, and full intensity usually set to 10 mW, small beam 1mm2,  
IsOptoSorted = 1;
[sData,filePath,fileName] = CreateSessionInfo(IsMoreStim,IsOptoSorted); %Log text file recorded during experiment will be opened and automatically looks for mouseinfo
%%% CONVERT LABVIEW RAW DATA INTO TDMS WHICH CAN BE USED BY MATLAB (Lick, PD, Water, Frame, Counter)
sData = loadTDMSdataNori(sData,filePath,fileName);
%%% calculate behavioral variables
%sData.daqdata.frameSignal(:,1) = 0;

BinNu = 80; % number of bins
%IsOptoSorted = 1; % was it an optical stimulated session? 0:no, 1:yes
SensitivityOpto = 5; % set 5 for sinus fT (also if stim goes to the next lap, until reward), set to 100 to short stimulations - then very sensitive
sData = CalcBehav(sData,BinNu,IsOptoSorted,SensitivityOpto);   % sData.behavior.opto.IsOptoSession = 1;
sData = OptoTrialSortingMoreStim(sData); % if more than one optical stimulation was used during session

%%% Comparing First - Second half - only behavior
%IsOptoSorted = 1; % if all trials was laser on or off
LapsTested = 20;
sData = BehavFirstSecondHalf(sData,LapsTested,IsOptoSorted); 

%%% ANALYSING INDIVIDUAL CALCIUM IMAGING SESSION:
% after extracting data from ROI manager you will have to open:
% _signals.mat_cia_deconvolved , _signals.mat_cia_denoised , _signals.mat_dff , _signals.mat_cia_spikes
% Use data extracted from RoiManager: open dff data, sData
VelMin = 0.1; %SET!!!! interneurons = 0 , PC= 0.1
pmtGainGreen = 10;
waveLength = 950;
laserPower = 125;              % mW
fovCoord = [0 0];          % AP, ML coordinated from center of window (0 -500 means imaging at -2.2 mm left hemisphere, ML +right, - left hem) left right changed 2019.04.03.
IsDeconv = 1; % 2: do the deconvolution, 1: use ROI manager deconvolution, 0 : do not use deconvolution
gol = 5; % gol#0 original reward, gol#2 50 cm forward, between hot glues, gol#3 50 cm before original after Velcro, gol = 10 no cues
sData = calcCaDataNori(sData,VelMin,laserPower,waveLength,fovCoord,pmtGainGreen,IsDeconv,gol); %calcCaDataNori %calcCaDataNoriShortVIP
clearvars IsOpto BinNu VelMin  laserPower waveLength fovCoord
close all
%Sorting Optodata and plotting heatmap and mean act for light-on/off/after trials:
%type = 1; %use dff, Slow Removed data
%type = 2; %use dff, Slow Removed, filtered data 
type = 0; %use dff unfiltered data
IsOptoSorted = 1; % if all the trials are the same =0, if there are light-on and off trials =1, if there are more opto protocol types= 2
sData = CaDataOptoSortingdFFsData(sData,type,IsOptoSorted);  
% if more than one optical stimulation was used during session
FigVisible = 'off';
sData  = CaDataOptoSortingdFFsDataMoreProts(sData,type,FigVisible); %type=0 dff, type=1 deconv data, type=2 spike rate 
sData  = CaDataAlignDeconvToStim(sData,FigVisible); 
sData  = CaDataAlignDeconvToLicks(sData,FigVisible); 
sData  = CaDataAlignDeconvToReward(sData,FigVisible); 

sData = CaDataonlySR(sData);
type = 1;
%%% PLace cell anal
InOutRatio = 3; %3 , 2
ReliabiliyIndex = 0.34; % 0.34 , 0.2
type = 0; % 0 : use dff data ,    -1 : use dff_SR data
sData = placeCellMaosData(sData,type,InOutRatio,ReliabiliyIndex);

% VIP data, recording the same cell over days
%ROIsDiscard = sData.imdata.ROIsToDiscard; % ROIs = []; analyse all ROIs
sData = placeCellMaosDataROIs(sData,type,InOutRatio,ReliabiliyIndex,sData.imdata.ROIsToDiscard);
close all;
 
%%% opto sessions:
sData = placecellMaosDataOptoOn(sData,type,InOutRatio,ReliabiliyIndex);
%%% more protocols
Protocol = 3;
sData = placecellMaosDataOptoOnMoreProt(sData,type,InOutRatio,ReliabiliyIndex,Protocol);  

%%% compare the same place cells in different sessions
plotmeanROIActAllSessions2(sData1,IsOptoSession,type) %%% have to set sctipt!

%%% compare the same place cells in different types of protocol within one session (light-on, light-off)
type = 0; % 0 dff % type = 1; % deconv
GaussFilter = 5;
sData = placecellComparison2(sData,type,GaussFilter);

newPlaceFields(sData)
