%% DATA CONVENTION
% Below is a suggestion to how we organize data in the lab. The most
% important thing is that we agree on the naming and structure. I have
% added suggestions. Make changes if you disagree and add stuff you think
% should be obligatory to put into the data. I have put the MATLAB class in
% parenteses behind all variables (char).

% I have used a data structure as a base for organizing the data which I
% call sData (short for session data). % If you are new to MATLAB read 
% about structs here: https://se.mathworks.com/help/matlab/ref/struct.html.
% I have also used CAPITALIZATION which is a common code convention for
% constants (which raw data sort of is). 

% My vision is that we make this file such that new people to the lab can 
% be given this file and start using the "correct" naming right away.

%% FUNCTION INPUT/OUTPUT ARGS
% When writing functions for whatever you do, I suggest that we use sData
% as the only input (and output if it is conventient). From there you can 
% extract whatever you need inside the functions. This makes it a lot 
% easier if you add / remove stuff from code later.

% Example: 
% function sData = my_function(sData)
% 
% lick_signal = sData.daqdata.lick_sig; % lick signal written compressed.
% 
% end


%% Organizing data
sData = struct(); % Main struct for storing all data. 

% 4 fields contain all the main data that I can think of now. New 
% avenues of data in the lab can be added if needed in the future.
sData.metadata = struct();
sData.daqdata = struct();
sData.imdata = struct();
sData.ephysdata = struct();

%% METADATA
% Hardwritten by user?
sData.metadata.SURGERY_DATE = '07.08.18'; % (char)
sData.metadata.ANIMAL_BREED = 'Wildtype'; % (char)
sData.metadata.VIRUS = 'AAV1.Syn.GCaMP7f'; % (char)
sData.metadata.EXPERIMENT_NOTES = 'I injected muscimol 2 hours before this recording.'; % (char)
sData.metadata.IMAGING_DEPTH = [130]; % In micrometer. If the piezo is used, upper and bottom depth should be used. 130 refers to -130 from bregma. Ex: for piezo over 120 ym it will be [130, 250]. (double)

? sData.metadata.LASER_POWER_MEASURED = 25; % in mW as measured by the user. (double)
? sData.metadata.LASER_WAVELENGTH_NM = 925; % in nanometer as set on the laser. (double)

% Read in from labview TDMS (Ideally your LabView code saves this for you)
sData.metadata.MOUSE_NAME = '1207'; % (char)
sData.metadata.SESSION_NUM = 1; % Session number (char)
sData.metadata.BLOCK_NUM = 1; % Block number (char)
sData.metadata.DAQ_SAMPLING_RATE = 5000; % (double)
sData.metadata.REC_DATE = '21.08.18'; % Date of the recording (char)
sData.metadata.REC_TIME = '12:41'; % Time of day when recording started (char)

% Read in from the imaging .ini file (this is already present in the .ini file)
sData.metadata.IMAGING = true; % true / false depending on if it is an imaging session or just raw behavioral data. If this says "false" every field below can be ignored in functions. (logical)
sData.metadata.MICROSCOPE = 'OS2'; % (char)
sData.metadata.XPIXELS = 512; % (double) 
sData.metadata.YPIXELS = 256; % (double)
sData.metadata.FOV_SIZE = '103 x 103 ym'; % in micrometer (char)
sData.metadata.IMAGING_FS = 30.86; % (double)
sData.metadata.ZOOM = 1.3750; % (double)
sData.metadata.N_CHANNELS = 2; % (double)
sData.metadata.CHANNEL_NAMES = {'Ch1','Ch2'}; % (cell)
sData.metadata.CHANNEL_COLOR = {'Green','Red'}; % (cell)
sData.metadata.POCKEL_VALUE = 37; % (double)

sData.metadata.PIEZO = 'NO'; % This should say NO if the piezo is not used, or SAW / ZIG depending on which mode is used if piezo is used.  (char)
sData.metadata.N_PLANES = 1; % (double)
sData.metadata.BESSEL = false; % true / false depending on whether the bessel beam module is used. (logical)

? sData.metadata.PSF_Z = 10; % Point spread function in Z? Might be handy if this bessel is used (double)

%% DAQDATA
% Typical signals obtained from the DAQ
sData.daqdata.LICK_SIG = []; % 1 x M array with each column being a sample. (double)
sData.daqdata.WHEEL_DIODE_SIG = []; % 1 x M array with each column being a sample. (double)
sData.daqdata.WATER_VALVE_SIG = []; % 1 x M array with each column being a sample. (double)
sData.daqdata.LFP_SIG = []; % 1 x M array with each column being a sample. (double)
sData.daqdata.FRAME_SIG = []; % The raw frame signal obtained from the microscope. 1 x N array with each column being a sample. (double)

? sData.daqdata.FRAME_ONSET = []; % 1 x N array with 1 for every sample where a new frame starts, and 0 elsewhere. I.e.: sum(sData.daqdata.FRAME_ONSET) will be equal to the number of frames.

% Processed signals
% I think that deciding on processed variable names is a bit too much. We
% should focus on raw data. However, somethings might be useful for us
% since so many are using running wheel setups. ALSO: Processed data will
% not be written CAPITALIZED. 

% For running wheel experiments:
sData.daqdata.wheel_position = []; % Position on the wheel from 0 to length of wheel ( ex 157.5 cm).
sData.daqdata.wheel_position_ds = []; % Same as above but downsampled to imaging frame rate.
sData.daqdata.lap = []; % Current lap from start of recording

%% IMDATA
sData.imdata.ROI_SIG_RAW = []; % N x M matrix with N ROIs and M samples (doubles)
sData.imdata.ROI_SIG_NEUROPIL = []; % N x M matrix with N ROIs and M samples (doubles)
sData.imdata.ROI_METADATA = struct(); % A struct containing all metadata for the ROIs. This is exported from roimanager. (struct)

sData.imdata.roi_sig_dff = []; % N x M matrix with N ROIs and M samples (doubles)
sData.imdata.roi_sig_deconv = []; % N x M matrix with N ROIs and M samples (doubles)
sData.imdata.roi_sig_dff_np_subtracted = []; % N x M matrix with N ROIs and M samples (doubles)


%% EPHYSDATA
% Anna; This is your area.
