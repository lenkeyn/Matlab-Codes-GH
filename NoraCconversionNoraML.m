% NoraC conversion to NoraML LV data

sData.behavior_old = sData.behavior;
sData.daqdata.waterValve = sData.daqdata.waterValveSignal;
sData.daqdata.optoSignal = 0;
fileID = sData.sessionInfo.sessionID;
sData.sessionInfo.fileID = sData.sessionInfo.sessionID;
savePath = 'C:\MATLAB\SAVE';
sData.sessionInfo.savePath = strcat(savePath,'\',fileID);
%save(fullfile(sData.sessionInfo.savePath,strcat(fileID,'_sData.mat')),'sData');

TDMSDATA.dt = 1/sData.daqdata.meta.fs; 
TDMSDATA.Distance = sData.daqdata.wheelRotaryEncoderSignal/sData.daqdata.meta.wheelRotaryEncoderTicks * 3.14*50;
TDMSDATA.Distance
% Dist = encoderSignal / tick * (pi*50)';
sData.daqdata.frameSignalOrig = sData.daqdata.frameSignal;
TDMSDATA.Framesignal = sData.daqdata.frameSignal;
TDMSDATA.Framesignal(TDMSDATA.Framesignal > mean(TDMSDATA.Framesignal))=1;
sData.daqdata.frameSignal = TDMSDATA.Framesignal;
FrameStart = diff(TDMSDATA.Framesignal)==1;
FrameStartIndex = find(FrameStart);

sData.daqdata.frameIndex = FrameStartIndex;         % double    Array of same lenght as frames required, where each sample n is the index number of other daqdata samples at the onset of frame n.
sData.daqdata.distanceCm = TDMSDATA.Distance;



%RUN
nBins = 80; % number of bins
IsOpto = 0; % was it an optical stimulated session? 0:no, 1:yes
OptoStimLimitMs = 15000; DiscardCmBeginningCm = 10; OptoSensitivity = 5; 
sData = CalcBehav2(sData,nBins,IsOpto,OptoSensitivity,OptoStimLimitMs,DiscardCmBeginningCm);
% except performance at end , comment out line 398-400
VelMin = 0.1; IsDeconv = 0;
sData = NoraCconversionNoraML_CaData(sData,VelMin,IsDeconv);

datatype = 0;
sData = placeCellMaosData2(sData,datatype);