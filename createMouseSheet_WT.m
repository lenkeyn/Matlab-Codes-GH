
SavePath = 'C:\MATLAB\MOUSEINFO';
mouseInfo = struct();                     % Struct containing all information about the mouse

% MOUSEINFO
mouseInfo.name = '8010';                                  % char      (REQUIRED) Mouse number according to lab convention
mouseInfo.dateOfBirth = '2018.04.00';                     % char      (REQUIRED) yyyy.mm.dd
mouseInfo.strain = 'B6';               % char      (REQUIRED) Strain of animal used. Be specific. 
mouseInfo.transgene = '-';
mouseInfo.Ordered = '-';
mouseInfo.sex = 'male';                                 % char      (REQUIRED) Male or female
mouseInfo.SLCageNumber = 31335;                           % double    (Optional) Science linker cage number
mouseInfo.SLMouseNumber = 82159;                          % double    (Optional) Science Linker animal number
%mouseInfo.SLMouseLitter = 0;
mouseInfo.lightCycle = 'Light on 10pm-10am';              % char      (Optional) Light is on during the following time, ex: '10 am to 10 pm'.

mouseInfo.surgeryDate = '2018.06.27';                          % char      (REQUIRED) yyyy.mm.dd
mouseInfo.windowCoordinates = 'AP -2.0 mm, ML 0 mm';      % double	'AP -2.2 mm, ML 0 mm'    (REQUIRED?) [x, y] (mm) distance from bregma to center of window.
mouseInfo.surgeryDoneBy = 'Nora';
mouseInfo.surgeryProtocol = 'virus injection + window implantation';                 %'Virus injection + window implantation'  char      (Optional) should refer to a documented protocol
mouseInfo.windowType = 'double, drilled whole on left';                 % char      Halfmoon/fullmoon
mouseInfo.injectedVirusN1 = 'AAV1-syn-GCamp7f';       % char      AAV1-CAG-flex-GCamp6s
mouseInfo.injectedVirusN1Location = 'Left and Right RSC';            %  'Left RSC'
mouseInfo.injectedVirusN1NanoLPerSite = 30;                % double
mouseInfo.injectedVirusN1NumberOfSites = 2;                % double
mouseInfo.injectedVirusN1Depth = 300;                      % double (micron)
mouseInfo.injectedVirusN2 = 'AAV8-syn-hM4D(Gi)-mCherry';       % char      AAV1-CAG-flex-GCamp6s
mouseInfo.injectedVirusN2Location = 'Left RSC';            %  'Left RSC'
mouseInfo.injectedVirusN2NanoLPerSite = 300;                % double
mouseInfo.injectedVirusN2NumberOfSites = 1;                % double
mouseInfo.injectedVirusN2Depth = 300;  
mouseInfo.RecordedHemisphere = 'right';                     % char      (Optional) relevant for injections/recordings
mouseInfo.injectionCoordinates = 'AP -2.0, ML -300';                      % double    (Optional) nInj x 2 (um) from bregma
mouseInfo.notes = 'NL-8, GCamp on both sides. 5 days control exp, then Muscimol or CTX Buffer injections to left RSC, imaging on right side, goal oriented spatial task. Skull was labelled for subicular injections.';

% save 
save(fullfile(SavePath,strcat('mouseInfo-',mouseInfo.name,'.mat')),'mouseInfo');
