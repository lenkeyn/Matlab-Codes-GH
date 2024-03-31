function = createsData

sData = struct();                               % Initialize the struct

sData.mouseInfo = struct();                     % Struct containing all information about the mouse
sData.sessionInfo = struct();                   % Struct containing information about the session 
sData.imdata = struct();                        % Imaging data
sData.daqdata = struct();                       % Signals obtained from the control computer DAQ
sData.ephysdata = struct();                     % Ephys signals
sData.behavior = struct();                      % Behavioral variables such as for example lick rate
sData.trials = struct();                        % Useful if the experiment contains trials.
sData.stats = struct();                         % Useful to add analysis variables and statistics


%% MOUSEINFO
sData.mouseInfo.name                            % char      (REQUIRED) Mouse number according to lab convention
sData.mouseInfo.dateOfBirth                     % char      (REQUIRED) yyyy.mm.dd
sData.mouseInfo.strain                          % char      (REQUIRED) Strain of animal used. Be specific. 
sData.mouseInfo.sex                             % char      (REQUIRED) Male or female
sData.mouseInfo.surgeryDate                     % char      (REQUIRED) yyyy.mm.dd
sData.mouseInfo.windowCoordinates               % double	(REQUIRED?) [x, y] (mm) distance from bregma to center of window.
sData.mouseInfo.lightCycle                      % char      (Optional) Light is on during the following time, ex: '10 am to 10 pm'.

sData.mouseInfo.SLCageNumber                    % double    (Optional) Science linker cage number
sData.mouseInfo.SLMouseNumber                   % double    (Optional) Science Linker animal number
sData.mouseInfo.surgeryProtocol                 % char      (Optional) should refer to a documented protocol
sData.mouseInfo.injectedVirus                   % char      (Optional) Relevant for injections
sData.mouseInfo.hemisphere                      % char      (Optional) relevant for injections/recordings
sData.mouseInfo.injectionCoordinates            % double    (Optional) nInj x 2 (um) from bregma


%% SESSIONINFO
sData.sessionInfo.sessionID                     % char      (REQUIRED) A unique ID for the session. It takes the form : mouseName-YYYYMMDD-SS (where SS is the session number).
sData.sessionInfo.date                          % char      (REQUIRED) yyyy.mm.dd
sData.sessionInfo.sessionNumber                 % double    (REQUIRED) The session number part of your sessionID.
sData.sessionInfo.sessionStartTime              % char      (REQUIRED) hh:mm:ss
sData.sessionInfo.sessionStopTime               % char      (REQUIRED) hh:mm:ss
sData.sessionInfo.recordedData                  % cell      (REQUIRED): Either {'2P','LFP','Patch','PupilVideo','SecurityCamera'}. These are only peripheral recordings for your DAQ, i.e. you do not need to add running wheel etc. 

sData.sessionInfo.mouseWeight                   % double    (REQUIRED, for water deprivation experiments)
sData.sessionInfo.mouseOriginalWeight           % double    (REQUIRED, for water deprivation experiments)
sData.sessionInfo.wasSessionAborted             % logic     (Optional) Was the session aborted before you wanted it to?
sData.sessionInfo.recDayNumber                  % double    (Optional) 1,2,3,...,N 
sData.sessionInfo.experimentName                % char      (Optional) name of experiment
sData.sessionInfo.protocol                      % char      (Optional) Tells where a documented protocol of your experiments are located. Fex Google docs etc.


%% DAQDATA
sData.daqdata.lickSignal                        % double    Lick signal
sData.daqdata.wheelRotaryEncoderSignal          % double    Rotary encoder signal used for running wheels
sData.daqdata.wheelDiode                        % double    Photodiode signal used to signify absolute position of wheel
sData.daqdata.waterValve                        % double    Water valve signal showing when the water valve is open
sData.daqdata.frameSignal                       % double    Frame signal recorded from the microscope
sData.daqdata.optoSignal                        % double    Optogenetics signal which is the voltage also sent to the laser diode 
sData.daqdata.frameIndex                        % double    Array of same lenght as frames required, where each sample n is the index number of other daqdata samples at the onset of frame n.


% DAQ Metadata
sData.daqdata.meta.fs                           % double    (REQUIRED) Sampling frequency of the data acquisition system
sData.daqdata.meta.wheelRotaryEncoderTicks      % double    (REQUIRED, if rotary encoder wheel is used) Number of ticks on the rotary encoder wheel used


%% IMDATA
sData.imdata.roiSignals = struct();             % struct    roiSignals is a struct array containing all signals as subfields. The channel color is specified as the array element number.
sData.imdata.roiSignals(2).ch                   % char      (REQUIRED) Use roiSignals(x).ch to signify the color of the channel. This is different between OsloScope I and II. Specify channel number in () and add .ch as 'red' or 'green'.

sData.imdata.roiSignals(2).roif                 % single    ROI fluroescence
sData.imdata.roiSignals(2).npilf                % single    Neuropil fluorescence
sData.imdata.roiSignals(2).dff                  % single    Delta F/F0

sData.imdata.roiSignals(2).deconv               % single    Deconvolved signal
sData.imdata.roiSignals(2).spikes               % single    Estimate spikes
sData.imdata.roiSignals(2).roifSubtractedNpil   % single    ROI fluorescence after subtracting neuropil
sData.imdata.roiSignals(2).dffSubtractedNpil    % single    Delta F/F0 after subtracting neuropil

sData.imdata.time                               % single    (Optional) If you want a time vector. This vector contains time stamps for each sample in roiSignals, starting at zero.


% Imaging Metadata (The following metadata can be import from the SciScan .ini file using the function getSciScanMetaData().)
sData.imdata.meta.microscope                    % char      (REQUIRED) Name of microscope
sData.imdata.meta.xpixels                       % double    (REQUIRED) Pixels in x dim
sData.imdata.meta.ypixels                       % double    (REQUIRED) Pixels in y dim
sData.imdata.meta.dt                            % double    (REQUIRED) Time difference between frames
sData.imdata.meta.fps                           % double    (REQUIRED) Frames per second (average)
sData.imdata.meta.zoomFactor                    % double    (REQUIRED) Zoom factor used for recording
sData.imdata.meta.zPosition                     % double    (REQUIRED) Relative depth of imaging in Z dim. The experimenter zeros the relative depth to dura.
sData.imdata.meta.umPerPxX                      % double    (REQUIRED) Micrometers per pixel in x dim
sData.imdata.meta.umPerPxY                      % double    (REQUIRED) Micrometers per pixel in y dim
sData.imdata.meta.nChannels                     % double    (REQUIRED) Number of recorded channels
sData.imdata.meta.channelNumbers                % double    (REQUIRED) Channels present, ex: [1,2].
sData.imdata.meta.channelNames                  % cell      (REQUIRED) Channel names, ex: {'Ch1','Ch2'}.
sData.imdata.meta.channelColor                  % cell      (REQUIRED) Channel colors, ex: {'Red','Green'}.
sData.imdata.meta.nFrames                       % double    (REQUIRED) Number of frames aquired
sData.imdata.meta.pockelSetting                 % double    (REQUIRED) Pockel setting in SciScan in percent
sData.imdata.meta.pmtGain                       % double    (REQUIRED) The gain set to PMT(s) during recording. (If two channels, [10, 20])
sData.imdata.meta.piezoActive                   % logical   (REQUIRED) true or false
sData.imdata.meta.piezoMode                     % char      (REQUIRED) Description of the mode used on the piezo, either SAW or ZIG
sData.imdata.meta.piezoNumberOfPlanes           % double    (REQUIRED) Number of planes used by the piezo.
sData.imdata.meta.piezoImagingRateHz            % double    (REQUIRED) The volume frames per second rate (basically overall imaging rate divided by number of piezo planes).
sData.imdata.meta.piezoVolumeDepth              % double    (REQUIRED) The number of micrometer the volume imaging spans.


% The below has to be handwritten by experimenter
sData.imdata.meta.laserPower                    % double    (REQUIRED) Laser power in mw
sData.imdata.meta.waveLength                    % double    (REQUIRED) Wavelength used on the laser
sData.imdata.meta.fovCoordinates                % double    (REQUIRED) [x, y, theta] micrometer between center of FOV and center of window. Theta is the rotation in degree.


% roiArray
sData.imdata.roiArray                           % struct    (REQUIRED) Struct array of same length as the number of ROIs. This is generated by roimanager and contains all metadata about the drawn ROIs, such as position in FOV, size etc. 


%% EPHYSDATA
sData.ephysdata.lfp = [];


% Meta
sData.ephysdata.meta.samplingRate                       % double    (REQUIRED) Sampling rate of ephys data
sData.ephysdata.meta.filterRange = [1000,100000];       % double    (REQUIRED) Upper and lower 
sData.ephysdata.meta.amplifier                          % char      (REQUIRED) Name and type of amplifier
sData.ephysdata.meta.amplifierGain                      % double    (REQUIRED) Gain on amplifier used
sData.ephysdata.meta.coordinateElectrodeTip             % double    (REQUIRED) [a.p, m.l., depth]
sData.ephysdata.meta.coordinateReferenceElectrodeTip    % double    (REQUIRED) [a.p, m.l., depth]


%% BEHAVIOR
% The naming of fields in behavior are as of yet up to each and every one
% to decide, but they should be perfectly explainable by reading the
% variable name. However, we advice you to use the following as
% abbreviations if needed:
% 
% Ds
%   All fields have a potential *Ds ending, which is ment for behavior
%   downsampled to the frames acquired from the microscope or potentially
%   some other signal, like LFP, tetrode recording, animal tracking etc.
% Binned
%   Binned may refer to spatial or temporal binning. In the case of a running
%   wheel, binning refers to spatial binning on the wheel with binSize being
%   in cm.
%
% Examples follow:

sData.behavior.wheelPos             % double    Absolute position on the wheel
sData.behavior.wheelPosDs           % double    Absolute position on the wheel (downsampled to the imaging FPS)
sData.behavior.wheelPosDsBinned     % double    Absolute position on the wheel (downsampled to the imaging FPS, and binned to some bin size stored in behavior.meta.binSize)

sData.behavior.wheelLap             % double    Absolute lap number from the start of the recording.

sData.behavior.runSpeed             % double    Running speed of the animal
sData.behavior.runSpeedDs           % double    Running speed of the animal (downsampled to the imaging FPS)

% Meta
sData.behavior.meta.binSize         % double    Bin size when binning data. For instance, this can be the spatial bin size used for binning linear track data.
sData.behavior.meta.nWheelBins      % double    Number of bins a running wheel has been divided into (useful for place field recordings)


%% TRIALS
% The fields of trials are up to everyone to decide. But make sure that the
% variable names are logical and easy to understand.
% sData.trials


%% STATS
% The field names used in stats are up to each individual to decide. But it
% is a conventient way of storing variables while you do analysis. 
% sData.stats