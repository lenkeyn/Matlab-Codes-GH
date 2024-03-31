
SavePath = 'C:\MATLAB\MOUSEINFO';
mouseInfo = struct();                     % Struct containing all information about the mouse

% MOUSEINFO
mouseInfo.name = '8025';                                  % char      (REQUIRED) Mouse number according to lab convention
mouseInfo.dateOfBirth = '2018.12.16';                     % char      (REQUIRED) yyyy.mm.dd
mouseInfo.strain = 'C57BL/6 x 129S4/SvJae';               % char      (REQUIRED) Strain of animal used. Be specific. 
mouseInfo.transgene = 'VIP-IRES-Cre';
mouseInfo.Ordered = 'https://www.jax.org/strain/010908';
mouseInfo.sex = 'female';                                 % char      (REQUIRED) Male or female
mouseInfo.SLCageNumber = 35870;                           % double    (Optional) Science linker cage number
mouseInfo.SLMouseNumber = 95065;                          % double    (Optional) Science Linker animal number
mouseInfo.SLMouseLitter = 12769;
mouseInfo.lightCycle = 'Light on 10pm-10am';              % char      (Optional) Light is on during the following time, ex: '10 am to 10 pm'.

mouseInfo.surgeryDate = '2019.01.24.';                    % char      (REQUIRED) yyyy.mm.dd
mouseInfo.windowCoordinates = 'AP -2.2 mm, ML 0 mm';      % double	(REQUIRED?) [x, y] (mm) distance from bregma to center of window.
mouseInfo.surgeryDoneBy = 'Koen';
mouseInfo.surgeryProtocol = 'Virus injection + window implantation';                 % char      (Optional) should refer to a documented protocol
mouseInfo.windowType = 'Halfmoon on left'; 
mouseInfo.injectedVirusN1 = 'AAV1-CAG-flex-GCamp6s';       % char      (Optional) Relevant for injections
mouseInfo.injectedVirusN1Location = 'Left RSC'; 
mouseInfo.injectedVirusN1NanoLPerSite = 15;                % double
mouseInfo.injectedVirusN1NumberOfSites = 5;                % double
mouseInfo.injectedVirusN1Depth = 200;                      % double
mouseInfo.RecordedHemisphere = 'Left';                     % char      (Optional) relevant for injections/recordings
%mouseInfo.injectionCoordinates = '';                      % double    (Optional) nInj x 2 (um) from bregma

% save 
save(fullfile(SavePath,strcat('mouseInfo-',mouseInfo.name,'.mat')),'mouseInfo');
