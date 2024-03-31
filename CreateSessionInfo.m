function [sData,filePath,fileName,stimProtocols] = CreateSessionInfo(IsMoreStim,IsOptoSorted)

%%% set SavePath
savePath = 'C:\MATLAB\SAVE';

%%% READ LAB BOOK DATA FTOM TXT FILE 

msgbox('Choose TDMSdata file');
[fileName,filePath,~] = uigetfile('*.tdms','','C:\' );

parts = strsplit(fileName,'.'); % Clear".tdms" from Filename.
sessionID = parts{1};

% Read lab book info to labBook 
if isfile([filePath sessionID '.txt'])
labBook = fileread([filePath sessionID '.txt']);
else
    msgbox('Lab book file containing meta data is required in .txt format with identical file name as the .tdms file.','Lab book file is not found.')
end

splitText = strsplit(labBook);
% Set indexes to extract data from splitText cell array
    weightIndex = 8;
    origWeightIndex = 15;
    weightPercIndex = 5;
    lapCompletedIndex = 55;

index = find(strcmp(splitText, 'Session'));
start = splitText(index(1)+2);
stop = splitText(index(1)+5);

index = find(strcmp(splitText, 'Weight:'));


%%% SESSIONINFO
sData.sessionInfo.sessionID = sessionID;
sData.sessionInfo.date = splitText{1};
sData.sessionInfo.sessionNumber = str2double(sessionID(16:17));
sData.sessionInfo.sessionStartTime = start{1};
sData.sessionInfo.sessionStopTime = stop{1};
sData.sessionInfo.lapCompleted = str2double(splitText{lapCompletedIndex});

% UPDATE LATER if there is imaging data
if sData.sessionInfo.sessionNumber < 1
   sData.sessionInfo.recordedData = {' '};  
else
   sData.sessionInfo.recordedData = {'2P'};
end

sData.sessionInfo.mouseWeight = str2double(splitText{weightIndex});
sData.sessionInfo.mouseOriginalWeight = str2double(splitText{origWeightIndex});
sData.sessionInfo.mouseWeightPercent = str2double(splitText{weightPercIndex});
sData.sessionInfo.labBook = labBook; 

% Test if the correct values are extracted
if isnan(sData.sessionInfo.mouseWeight) + isnan(sData.sessionInfo.mouseOriginalWeight) + isnan(sData.sessionInfo.mouseWeightPercent) > 0
    msgbox('Modify indexes to extract correct values from splitText cell array.','Incorrect indexes!')
    return
end

%%% MOUSEINFO
mouseFolder = 'C:\MATLAB\MOUSEINFO';

if isfile([mouseFolder '\mouseInfo-' sessionID(2:5) '.mat'])
    load([mouseFolder '\mouseInfo-' sessionID(2:5)]); % fileName might change
    sData.mouseInfo = mouseInfo;
    clear('mouseInfo');
else
    msgbox(['1) Make sure if "mouseFolder" variable is defined in the script and refers to the path where the mouseinfo files are stored.',char(10),char(10),... 
        '2) Make sure if the mouse info file for this mouse is filled in manually and saved according to the mouse naming standards.'],'Mouseinfo data is not found.');
end



%%% optical stimulation protocol
if IsOptoSorted == 1 %IsMoreStim == 1
    
    % Extract max light intensity
    if size(strsplit(labBook,'Photo stimulation:'),2) > 1
        splitText = strsplit(labBook,'Photo stimulation:');
        splitText = strsplit(splitText{2},'*****');
        splitText = strsplit(splitText{1},'Max output: ');
        splitText = strsplit(splitText{2},'light intensity: ');
        splitText = strsplit(splitText{2},'\r\n');
        text = splitText{1}; 
        if strfind(text,'mW/mm2') > 1
            maxStimInt_mWpermm2 = str2double(text(1:strfind(text,'mW/mm2')-1));

        elseif strfind(text,'mW') > 1
            maxStimInt_mW = str2double(text(1:strfind(text,'mW')-1));

            % get illuminated surface area
            splitText = strsplit(labBook,'surface area: ');
            if size(splitText,2) > 1
                splitText = strsplit(splitText{2},'mm2');
                beamSize = str2double(splitText{1});
            else        
                beamSize = inputdlg('Beam size of the stimulus laser at the working distance of the objective (mm2): ','Illuminated surface area was not find in the labbok. Enter the value manually!');
                beamSize = str2double(beamSize{1});
                maxStimInt_mWpermm2 = maxStimInt_mW / beamSize;
            end

        else
            msgbox('The stimulus intensity value is not given properly in the lab book file. The units should be either "mW" or "mW/mm". The intensity should be a numerical value.','Error!');
        end
    end

% Extract optical stimulation protocol
    splitText = strsplit(labBook,'Optical stimulation protocol:');
    splitTextM = strsplit(splitText{1},'Primary masking light:');
    splitText = strsplit(splitText{2},'Masking protocol (secondary masking light):');
    splitText = strsplit(splitText{1},'\r\n');
    stimProtocols(1:numel(splitText)-2) = struct();
    for i = 1:1:numel(splitText)-2
       % trialTypeIndicator
       splitText2 = strsplit(splitText{i+1},': ');
       stimProtocols(i).trialTypeIndicator = str2double(splitText2{1});
       stimProtocols(i).trialType = i-1;
       % proportion
       splitText2 = strsplit(splitText{i+1});
       splitText2 = strsplit(splitText2{2},'/');
       stimProtocols(i).proportion = str2double(splitText2{1}) / str2double(splitText2{2});
       % full protocol
       splitText2 = strsplit(splitText{i+1},['/' splitText2{2} ' ']);
       stimProtocols(i).protocol = splitText2{2};
       if ~isequal(stimProtocols(i).protocol,'none')
           % waveform
           splitText2 = strsplit(stimProtocols(i).protocol,' from ');
           stimProtocols(i).waveform = splitText2{1};
           % from
           splitText2 = strsplit(splitText2{2},' to ');
           if isequal('rew',splitText2{1})
               stimProtocols(i).from = 'reward';
           else
               stimProtocols(i).from = str2double(splitText2{1});
           end
           % to
           splitText2 = strsplit(splitText2{2},' int. ');
           if isequal('rew',splitText2{1})
               stimProtocols(i).to = 'reward';
           else
               stimProtocols(i).to = str2double(splitText2{1});
           end
           % intensity
           splitText2 = strsplit(splitText2{2},' %');
           stimProtocols(i).voltageRangePerc = str2num(splitText2{1});
           %stimProtocols(i).intensity_mWpermm2 = stimProtocols(i).intensity_perc/100*maxStimInt_mWpermm2;
           % intensity in mW/mm2. Fill out manually later if you want
           %splitText2 = strsplit(splitText2{2},' %');
           %stimProtocols(i).intensity_mWmm2 = str2num(splitText2{1});
           % masking light (primary)
           stimProtocols(i).maskingLight = splitTextM{2}(1:end-4);
       end
    end
    stimProtocols(1).maskingLight = splitTextM{2}(1:end-4); % masking light for control

    sData.stimProtocols = stimProtocols;
    clear('stimProtocols');
end
    
%%% SAVING
parts = strsplit(sessionID,'-'); % short version of sessionID
fileID = strcat(parts{1},'-',parts{2},'-',parts{3});
sData.sessionInfo.fileID = fileID;
mkdir(savePath,fileID);
sData.sessionInfo.savePath = strcat(savePath,'\',fileID);
save(fullfile(sData.sessionInfo.savePath,strcat(fileID,'_sData.mat')),'sData');

end