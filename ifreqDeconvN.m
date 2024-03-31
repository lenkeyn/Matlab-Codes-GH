function sData = ifreqDeconvN(sData)
% calculates event weighted frequency of events that are separated by 0s
% 
%signal = ciaDeconv(roi,:);

% for already extracted deconvolved signal
fps = sData.imdata.meta.fps;
signal = sData.imdata.roiSignals(2).deconv;
nROIs = size(signal,1);
nSamples = size(signal,2);

iFrequency = nan(size(signal));

for roi = 1:1:nROIs
    ROIsignal = signal(roi,1:nSamples);
    eventPos = find(ROIsignal ~= 0);
    events = ROIsignal(eventPos);                          % first datapoint is excluded
    freques = fps./diff(eventPos);
    freques(size(freques,2)+1) = 0;    % add 0 as the last datapoint because frequency can only be calculated berween two events 
    ROIiFrequency = nan(1,nSamples);
    ROIiFrequency(ROIsignal ~= 0) = freques.*events;          % fill data in the correct positions 
    ROIiFrequency(1) = 0;                                  % add 0 as the first value for the fillmissing function
    ROIiFrequency = fillmissing(ROIiFrequency,'previous'); 
    % divide data with the smalest "unitary" deconvolved value
    if max(signal(1:nSamples)) > 0
        ROIiFrequency(1:nSamples) = ROIiFrequency(1:nSamples) / min(signal((signal(1:nSamples) > 0))); 
    end    
    iFrequency(roi,1:nSamples) = ROIiFrequency(1:nSamples);
end

sData.imdata.roiSignals(2).iFrequency = single(iFrequency);

end