function sDataFiles = loadDataMate(method,filePath)
%%% Load multiple sData files for analysis

% If method = 'light' clears daqdata and ROI signals for faster operation
% and saving RAM memory.

if ~exist('method')
    method = '';
end

%if nargin < 2
% select single or multiple files to analyze 
%[fileNames,filePath,~] = uigetfile('*.mat','','C:\Users\Mate Neubrandt\Documents\RECORDINGS','MultiSelect','on' );
%end

if ~iscell(filePath)
    filePath = {filePath};
end
fileNumber = length(filePath);
sDataFiles = cell(1,fileNumber);

for f = 1:1:fileNumber
    fileName = filePath{f};
    load(fileName);
    % insert analysis function here:
    if strcmp(method,'light')
        sData = rmfield(sData,'daqdata');
    end
    if isfield(sData,'imdata') 
        sData.imdata = rmfield(sData.imdata,'roiSignals');
    end
    sDataFiles{f} = sData;
    clear('sData');    
end

end