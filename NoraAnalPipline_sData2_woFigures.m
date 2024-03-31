%%% NORI ANALYSIS PIPLINE - SDATA ERA

%%%%% ANALYSIS OF IMAGED CA-TRANSIENTS AND BEHAVIOR
%{
%%% Create Session Info, DAQData Info, Behavioral Info:
IsMoreStim = 0; % were there different opto-stim protocols?
RedLaserMax = 10; % 100% laser stim mW/mm2. Big beam: 10-11 mm2, and full intensity usually set to 10 mW, small beam 1mm2, full intensity was set to 10 mW/mm2 
IsOptoSorted = 1;
[sData,filePath,fileName] = CreateSessionInfo(IsMoreStim,IsOptoSorted); %Log text file recorded during experiment will be opened and automatically looks for mouseinfo
%%% CONVERT LABVIEW RAW DATA INTO TDMS WHICH CAN BE USED BY MATLAB (Lick, PD, Water, Frame, Counter)
sData = loadTDMSdataNori(sData,filePath,fileName);
%%% calculate behavioral variables
%sData.daqdata.frameSignal(:,1) = 0;
%}

nBins = 80; % number of bins
IsOpto = 1; OptoSensitivity = 5; OptoStimLimitMs = 9000; DiscardCmBeginningCm = 10; 
sData = CalcBehav2(sData,nBins,IsOpto,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm); close all; %sData = CalcBehav(sData,BinNu,IsOptoSorted,SensitivityOpto);   % sData.behavior.opto.IsOptoSession = 1;
% more stim
%sData = OptoTrialSortingMoreStim2(sData); % if more than one optical stimulation was used during session. Failed trials labelled with -1, dark grey on plot.

%%% Comparing First - Second half - only behavior
LapsTested = 20; 
sData = BehavFirstSecondHalf(sData,LapsTested,IsOpto); close all;

%%% ANALYSING INDIVIDUAL CALCIUM IMAGING SESSION:
VelMin = 0.1; %SET!!!! interneurons = 0 , PC= 0.1
IsDeconv = 1; % 2: do the deconvolution, 1: use ROI manager deconvolution, 0 : do not use deconvolution
gol = 5; % gol#0 original reward, gol#2 50 cm forward, between hot glues, gol#3 50 cm before original after Velcro, gol = 10 no cues
sData = calcCaDataNori2_woFigures(sData,VelMin,IsDeconv,gol); %calcCaDataNori %calcCaDataNoriShortVIP
clearvars  BinNu VelMin  
close all
%Sorting Optodata and plotting heatmap and mean act for opto-on/off/after trials:
FigVisible = 'off';
sData = CaDataOptoSortingdFFsData2(sData,FigVisible); %sData  = CaDataOptoSortingDeconvData2(sData); % no smoothing %sData  = CaDataOptoSortingCumSpikesData2(sData); % no smoothing

%%% Place cells in optically stimulated sessions:
datatype = 0; % use dff data , datatype = 1 for deconvolved data
sData = placecellMaosDataOptoOn2(sData,datatype,FigVisible);
sData = LandmarkCellDetectionOptoOn(sData); % detects cells which respond to both cues (e.g. both hot glue spikes or both velcro), detects place cells with only one place field as well

%%% compare the same place cells in different types of protocol within one session (opto-on, opto-off)
type = 0; % 0 dff % type = 1; % deconv
GaussFilter = 5;
sData = placecellComparison3(sData,type,GaussFilter);

%%% new place field analysis
newPlaceFields(sData)


% gain modulation scripts
sData.gainModulationPCinAllProtocol = struct;
sData.gainModulationPCinAnyProtocol = struct;
DiscardCmBeginning = 20;
figGeneration = 1; % 1: on, 0 = off
sData = gainMod_binnedOptoOffOnRegressionPCinAllProtocols(sData,figGeneration,DiscardCmBeginning);
figGeneration = 1;
sData = gainMod_binnedOptoOffOnRegressionPCinAnyProtocols(sData,figGeneration,DiscardCmBeginning);
sData = gainModulationInOutRatioPCinAllProtocols(sData);
sData = gainModulationInOutRatioPCinAnyProtocols(sData);

DiscardCmBeginning = 20;
figGeneration = 1;
sData.gainModulationPCinCtr = struct;
sData.gainModulationPCinLandmarkAndOnePFPC = struct;
sData = gainMod_binnedOptoOffOnRegressionPCLandmarkAndOnePFPC(sData,figGeneration,DiscardCmBeginning);
sData = gainMod_binnedOptoOffOnRegressioninCtr(sData,figGeneration,DiscardCmBeginning);
sData = gainModulationInOutRatioPCinCtr(sData);
sData = gainModulationInOutRatioLandmarkAndOnePFPC(sData);

sumTable2 = sumImagingOpto2(sData);
sumTable3 = sumImagingOpto3(sData);

% Save file to same path where other files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
