function sData = loadTDMSdataNori(sData,filePath,fileName)

% OPEN TDMS DATA FILE WITH LICK-PD-WATERWALVE-FRAME DATA, LOAD DATA
% [fileName,filePath,~] = uigetfile('*.tdms');

channelNames = {'Distance','Framesignal','Photodiode','Valveopen','Licksignal','Photostimulussignal','Encodersignal','Stimulusprotocol'}; % processed data
%changed strange values to nothing, because Matlab cannot handle it
channelNames = cellfun( @(chName) regexprep(chName, '_', ''), channelNames, 'uni', false);
channelNames = cellfun( @(chName) regexprep(chName, '__', ''), channelNames, 'uni', false);
channelNames = cellfun( @(chName) regexprep(chName, '___', ''), channelNames, 'uni', false);

% % Convert tdms-file with daq data to mat-file
matFile = simpleConvertTDMS(fullfile(filePath,fileName));

% % Load matfile and create an array from datastruct fields.
S1 = load(matFile{1});
fieldNames = fields(S1); % S1 data file column names are put in an array

% Create a struct to hold the data.
TDMSDATA = struct;

% Go through each channel and add data to the data struct
for i = 1:numel(channelNames)
    varName = channelNames{i}; 
    fieldPos = find((contains(fieldNames, varName)),1); % find the postition of the ChannelNames in the fieldNames array, and put fieldPos array. (There will be 2 position, but only one will contain data)
    if ~isempty(fieldPos) % if filedPos is not empty = varNames were found in fieldNames array
        for j = 1:numel(fieldPos) % do as many times as varName can be found in fieldNames
            S1data = S1.(fieldNames{fieldPos(j)}).Data; % get the data from S1 proper column
            if isempty(S1data) % if it is empty (=1), discard and continue to find the data, jumps back to beginning of for cycle.
                continue
            else  % if it is not empty, put the data into struct
                % struct fields cannot start with a number...
                % if strcmp(varName, '2PFrames') % string-compare. Returns =1 if the two strings are identical
                %    varName = 'FrameSignal2P';
                % end
                TDMSDATA.(varName) = S1data;
            end
        end
    else
        %warning('Channel name ''%s'' was not found in tdms file ''%s''', varName, filename)
    end
end

% Add delta t and sampling rate to struct. Only works for TDMS files 
TDMSDATA.dt = 1/3000; % sampling frequency was set to 3 kHz from 2018.10.17., before it was 200 Hz
TDMSDATA.Distance(1:60)=0; % there is a crazy number

%TDMSDATA.Framesignal(1:end)=0;
if sum(TDMSDATA.Framesignal) == 0
    TDMSDATA.Fakeframesignal = TDMSDATA.Framesignal;
    Samples = 1:1:numel(TDMSDATA.Framesignal);
    TDMSDATA.Fakeframesignal(rem(Samples,97)==0) = 1;
end

TDMSDATA.filename = fileName;
TDMSDATA.filePath = filePath;

%%% calculate frameStart indices
FrameStart = diff(TDMSDATA.Framesignal)==1;
FrameStartIndex = find(FrameStart);

%%% generate sDATA.DAQDATA

sData.daqdata = struct();

sData.daqdata.lickSignal = TDMSDATA.Licksignal;     % double    Lick signal                       
sData.daqdata.wheelRotaryEncoderSignal = TDMSDATA.Encodersignal; % double    Rotary encoder signal used for running wheels
sData.daqdata.wheelDiode = TDMSDATA.Photodiode;     % double    Photodiode signal used to signify absolute position of wheel
sData.daqdata.waterValve = TDMSDATA.Valveopen;      % double    Water valve signal showing when the water valve is open
sData.daqdata.frameSignal = TDMSDATA.Framesignal;   % double    Frame signal recorded from the microscope
sData.daqdata.frameSignal(1:1000,1) = 0; % many times there is an artefact-framesignal in the beginning, when the imaging surely does not start yet
sData.daqdata.optoSignal = TDMSDATA.Photostimulussignal;   % double    Optogenetics signal which is the voltage also sent to the laser diode 
if any(ismember(fields(TDMSDATA),'Stimulusprotocol')) &&  sum(TDMSDATA.Stimulusprotocol) > 0
    sData.daqdata.optoStimProtocol = TDMSDATA.Stimulusprotocol; %double , if there were more opto-protocols
end
sData.daqdata.frameIndex = FrameStartIndex;         % double    Array of same lenght as frames required, where each sample n is the index number of other daqdata samples at the onset of frame n.
sData.daqdata.distanceCm = TDMSDATA.Distance;
% DAQ Metadata
sData.daqdata.meta.fs = 3000;                           % double    (REQUIRED) Sampling frequency (Hz) of the data acquisition system
sData.daqdata.meta.wheelRotaryEncoderTicks = 2000;      % double    (REQUIRED, if rotary encoder wheel is used) Number of ticks on the rotary encoder wheel used
sData.daqdata.notes = 'Encoder signal was converted to Distance(cm) during recording: Dist = encoderSignal / tick * (pi*50)';


% Save file to same path where LV files can be found 
save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end