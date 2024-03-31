function sData = LandmarkCellDetectionOptoOnDFFVIPpaper3(sData)

%%% CHECK FOR LANDMARK CELLS: 
% there are two set of landmarks on the wheel at certain distance: 
% velcro at 26 and 86 cm and hot glue spikes at 46 and 66 cm. All tactile cues span for 2 cms.
% The distance between HGS are 20 cm and velcros are 60 cm. I look for double peaks with these distances. 
% Since in the dF/F data the second peak is often masked by the first peak's decay, I check these double peaks in the dconvolved (unsmoothed) data

% I will check landmark properties in each opto protocol (off, on, after)

LandmarkArea = floor(18/sData.behavior.meta.binSize); % 20 cm (changed from 4 cm), search for landmark cells where the landmark was presented and a little addition
VelcroStart1 = ceil(26/sData.behavior.meta.binSize); % 
VelcroEnd1 = VelcroStart1 + LandmarkArea;
VelcroStart2 = floor(86/sData.behavior.meta.binSize);
VelcroEnd2 = VelcroStart2 + LandmarkArea;
HGSStart1 = ceil(46/sData.behavior.meta.binSize);
HGSEnd1 = HGSStart1 + LandmarkArea;
HGSStart2 = ceil(66/sData.behavior.meta.binSize);
HGSEnd2 = HGSStart2 + LandmarkArea; % delay in Ca transient due to previous Ca transient evoked by first HGS
DiffVelcro = floor(60/sData.behavior.meta.binSize); %30 bin
DiffHGS = floor(20/sData.behavior.meta.binSize); %10 bin
AllowedJitter = 1; % bins
MinPeakDetectionMainPeak = sData.imdata.roiStat.meanPeakDff/10; % spike detection treshold in dFF, it was 0.01 , 0.1 is too much
%MinPeakDetectionMainPeak = 0.07;
%MinPeakDetectionSecondPeak = sData.imdata.roiStat.meanPeakDff/20; % spike detection treshold in dFF, it was 0.0003

datasetDFF = sData.imdata.binned.RoidFF;
LandmarkDetectionDFF = struct;

LandmarkCells = struct;    
LandmarkCells.LandmarkCells = NaN(sData.imdata.nROIs,3); 
LandmarkCells.VelcroCells = NaN(sData.imdata.nROIs,3);
LandmarkCells.HGSCells = NaN(sData.imdata.nROIs,3);
LandmarkCells.Peak1PosBin = NaN(sData.imdata.nROIs,3);
LandmarkCells.Peak2PosBin = NaN(sData.imdata.nROIs,3);

for k = 1:1:3 % run in each protocol
    
    if k == 1
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoOff; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.OptoOffTrials;  
    elseif k == 2
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoOn; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.OptoOnTrials; 
    elseif k == 3    
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoAfter; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.AfterOptoTrials;   
    end 
    
    for i = 1:1:sData.imdata.nROIs % comparing peak positions (for a ROI) if the distance of peaks equals to distance between landmarks 
        if PCdff.Criteria123Passed(i) > 0 % if the cell was a place cell in that protocol check if it was a landmark cell
            DataTemp1 = [nanmean(datasetDFF{1,i}(Trials,:),1) nanmean(datasetDFF{1,i}(Trials,:),1)];  % concatentate the data twice to be able to use as circular dataset
            [pks,locs] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionMainPeak); % peaks and indices at which the peaks occur. Very local peaks are also detected
            PeakAmplitudeLocation = locs(find(pks == max(pks),1)); % search for the peak amplitude's location (which bin)
            PeakAmplitude = max(pks); 
            if isempty(PeakAmplitudeLocation)
                continue
            end
            % check if it can be a velcro or a HGS cell, then check if there is a peak at the other similar landmark (one bin jitter allowed for Velcro, 2 bin for HGS)
            if PeakAmplitudeLocation >= VelcroStart1 &&  PeakAmplitudeLocation <= VelcroEnd1 % Velcro1 is the main peak, the second peak should be at least this large
                MinPeakDetectionSecondPeak = PeakAmplitude/5;
                DiffVelcroArray = PeakAmplitudeLocation + DiffVelcro - AllowedJitter : PeakAmplitudeLocation + DiffVelcro + AllowedJitter;
                DiffVelcroArray(DiffVelcroArray < VelcroStart2) = NaN;
                DiffVelcroArray(DiffVelcroArray > VelcroEnd2) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if  any(ismember(locs2,DiffVelcroArray)>0)
                    LandmarkCells.LandmarkCells(i,k) = i; 
                    LandmarkCells.VelcroCells(i,k) = i;
                    LandmarkCells.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffVelcroArray)>0));
                end
            elseif PeakAmplitudeLocation >= VelcroStart2 &&  PeakAmplitudeLocation <= VelcroEnd2 % Velcro2 is the main peak
                MinPeakDetectionSecondPeak = PeakAmplitude/5;
                DiffVelcroArray = PeakAmplitudeLocation - DiffVelcro - AllowedJitter : PeakAmplitudeLocation - DiffVelcro + AllowedJitter;
                DiffVelcroArray(DiffVelcroArray < VelcroStart1) = NaN;
                DiffVelcroArray(DiffVelcroArray > VelcroEnd1) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffVelcroArray)>0) 
                    LandmarkCells.LandmarkCells(i,k) = i; 
                    LandmarkCells.VelcroCells(i,k) = i;
                    LandmarkCells.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffVelcroArray)>0));
                end
            elseif PeakAmplitudeLocation >= HGSStart1 &&  PeakAmplitudeLocation <= HGSEnd1 % HGS1 is the main peak
                MinPeakDetectionSecondPeak = PeakAmplitude/10; % because dFF tail of the first peak often still there
                DiffHGSArray = PeakAmplitudeLocation + DiffHGS - AllowedJitter : PeakAmplitudeLocation + DiffHGS + AllowedJitter + 1; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
                DiffHGSArray(DiffHGSArray < HGSStart2) = NaN;
                DiffHGSArray(DiffHGSArray > HGSEnd2) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffHGSArray)>0) 
                    LandmarkCells.LandmarkCells(i,k) = i; 
                    LandmarkCells.HGSCells(i,k) = i;
                    LandmarkCells.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffHGSArray)>0));
                end
            elseif PeakAmplitudeLocation >= HGSStart2 &&  PeakAmplitudeLocation <= HGSEnd2 % HGS2 is the main peak
                MinPeakDetectionSecondPeak = PeakAmplitude/5; % because dFF tail of the first peak often still there
                DiffHGSArray = PeakAmplitudeLocation - DiffHGS - AllowedJitter -1 : PeakAmplitudeLocation - DiffHGS + AllowedJitter; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
                DiffHGSArray(DiffHGSArray < HGSStart1) = NaN;
                DiffHGSArray(DiffHGSArray > HGSEnd1) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffHGSArray)>0)
                    LandmarkCells.LandmarkCells(i,k) = i; 
                    LandmarkCells.HGSCells(i,k) = i;
                    LandmarkCells.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffHGSArray)>0));
                end
            end
        end
    end     
end

for m = 1:1:3
    LandmarkCellListPre = LandmarkCells.LandmarkCells(1:end,m);
    VelcroCellsListPre =  LandmarkCells.VelcroCells(1:end,m);
    HGSCellsListPre = LandmarkCells.HGSCells(1:end,m);
    Peak1PosBinPre = LandmarkCells.Peak1PosBin(1:end,m);
    Peak2PosBinPre = LandmarkCells.Peak2PosBin(1:end,m);
    if m == 1
        LandmarkDetectionDFF.OptoOff.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetectionDFF.OptoOff.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetectionDFF.OptoOff.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetectionDFF.OptoOff.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetectionDFF.OptoOff.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    elseif m == 2
        LandmarkDetectionDFF.OptoOn.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetectionDFF.OptoOn.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetectionDFF.OptoOn.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetectionDFF.OptoOn.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetectionDFF.OptoOn.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    elseif m == 3
        LandmarkDetectionDFF.OptoAfter.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetectionDFF.OptoAfter.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetectionDFF.OptoAfter.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetectionDFF.OptoAfter.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetectionDFF.OptoAfter.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    end
end
       

LandmarkDetectionDFF.Params.LandmarkArea = LandmarkArea; 
LandmarkDetectionDFF.Params.VelcroStart1 = VelcroStart1;
LandmarkDetectionDFF.Params.VelcroEnd1 = VelcroEnd1;
LandmarkDetectionDFF.Params.VelcroStart2 = VelcroStart2;
LandmarkDetectionDFF.Params.VelcroEnd2 = VelcroEnd2;
LandmarkDetectionDFF.Params.HGSStart1 = HGSStart1;
LandmarkDetectionDFF.Params.HGSEnd1 = HGSEnd1;
LandmarkDetectionDFF.Params.HGSStart2 = HGSStart2;
LandmarkDetectionDFF.Params.HGSEnd2 = HGSEnd2;
LandmarkDetectionDFF.Params.DiffVelcro = DiffVelcro;
LandmarkDetectionDFF.Params.DiffHGS = DiffHGS;
LandmarkDetectionDFF.Params.AllowedJitter = AllowedJitter; 
LandmarkDetectionDFF.Params.MinPeakDetection = MinPeakDetectionSecondPeak;


LandmarkDetectionDFF.LandmarkCellinAllProt = intersect(intersect(LandmarkDetectionDFF.OptoOff.LandmarkCellsList,LandmarkDetectionDFF.OptoOn.LandmarkCellsList),LandmarkDetectionDFF.OptoAfter.LandmarkCellsList);%,MaoOpto.OptoAfter.PlaceCells);
LandmarkDetectionDFF.LandmarkCellinAnyProt = union(union(LandmarkDetectionDFF.OptoOff.LandmarkCellsList,LandmarkDetectionDFF.OptoOn.LandmarkCellsList),LandmarkDetectionDFF.OptoAfter.LandmarkCellsList);%,MaoOpto.OptoAfter.PlaceCells);
LandmarkDetectionDFF.LandmarkCellinOptoOfforOptoOn = union(LandmarkDetectionDFF.OptoOff.LandmarkCellsList,LandmarkDetectionDFF.OptoOn.LandmarkCellsList);

sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF = LandmarkDetectionDFF;
        

%%% detect place cells with only one place field
%datasetDFF = sData.imdata.binned.RoidFF; % use dFF data
    
PlaceCellsWithOnePF = struct;
PlaceCellsWithOnePF.PlaceCellsWithOnePF = NaN(sData.imdata.nROIs,3);
PlaceCellsWithOnePF.PeakPosdFFBin = NaN(sData.imdata.nROIs,3);
PlaceCellsWithOnePF.PeakPosDeconvBin = NaN(sData.imdata.nROIs,3);

for k = 1:1:3 % run in each protocol
    
    if k == 1
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoOff; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.OptoOffTrials;  
    elseif k == 2
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoOn; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.OptoOnTrials; 
    elseif k == 3    
        PCdff = sData.imdata.MaoPC_Opto_dff.OptoAfter; % struct containing which cells are place cells based on dF/F data
        Trials = sData.imdata.binned.AfterOptoTrials;   
    end 
    
    for i = 1:1:sData.imdata.nROIs % search for place cells with one peak 
        if PCdff.Criteria123Passed(i) > 0 && ~any(sData.imdata.MaoPC_Opto_dff.LandmarkCellsDFF.LandmarkCellinOptoOfforOptoOn==i)  % if the cell is a place cell but not a landmark cell
            %DataTemp = smoothdata(nanmean(dataset{1,i},1),'Gaussian',3);  % use little bit of smoothing on deconvolved data, because sometimes larger transients cause double peaks in deconvolution
            %[pks,locs] = findpeaks(DataTemp,'MinPeakProminence',MinPeakDetection); % peaks and indices at which the peaks occur.
            DataTemp2 = [nanmean(datasetDFF{1,i}(Trials,:),1) nanmean(datasetDFF{1,i}(Trials,:),1)]; 
            DataTempMax = max(DataTemp2);
            MinPeakDetectionSecondPeak = DataTempMax/5;
            [pks,locs] = findpeaks(DataTemp2,'MinPeakProminence',MinPeakDetectionSecondPeak); % peaks and indices at which the peaks occur.
            PeakAmplitude = max(pks);
            PeakAmplitudeLocation = locs(find(pks == max(pks),1)); % search for the peak amplitude's location (which bin)
            %if length(PeakAmplitudeLocation) < 1
            %    continue
            %end
            pksMod = pks;
            pksMod(pksMod == PeakAmplitude) = 0;
            if any((pksMod)> PeakAmplitude/3) % if there is another big peak in the dFF data (1/3 of biggest), it is not considered as a single peak place cell
                continue
            else
                PlaceCellsWithOnePF.PlaceCellsWithOnePF(i,k) = i;
                PlaceCellsWithOnePF.PeakPosdFFBin(i,k) = PeakAmplitudeLocation;
                if PeakAmplitudeLocation > sData.behavior.meta.nBins
                    PlaceCellsWithOnePF.PeakPosdFFBin(i,k) = PeakAmplitudeLocation - sData.behavior.meta.nBins;
                end
                % find the peak location in the deconvolved signal
                DataTemp3 = nanmean(datasetDFF{1,i},1);
                DataTempMax = max(DataTemp3);
                MinPeakDetectionSecondPeak = DataTempMax/4;
                [~,locs] = findpeaks(DataTemp3,'MinPeakProminence',MinPeakDetectionSecondPeak); % peaks and indices at which the peaks occur.
                if isempty(locs) 
                    locs = find(DataTemp3 == DataTempMax);
                end
                locs2 = locs(locs <= PeakAmplitudeLocation);
                if isempty(locs2)
                   locs2 = locs;
                   PeakAmplitudeLocation = PeakAmplitudeLocation + sData.behavior.meta.nBins; % if the dff peak is in the beginning of trial and deconv peak is at the end of previous trial
                end
                locs3 = PeakAmplitudeLocation - locs2;
                PlaceCellsWithOnePF.PeakPosDeconvBin(i,k) = locs2(locs3 == min(locs3));
            end
        end
    end
    %PlaceCellsWithOnePF.PlaceCellsWithOnePFList = PlaceCellsWithOnePF.PlaceCellsWithOnePF(~isnan(PlaceCellsWithOnePF.PlaceCellsWithOnePF));
end
% write data into struct
PlaceCellsWithOnePFOff = PlaceCellsWithOnePF.PlaceCellsWithOnePF(:,1);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOff = PlaceCellsWithOnePFOff(~isnan(PlaceCellsWithOnePFOff));
PlaceCellsWithOnePFOn = PlaceCellsWithOnePF.PlaceCellsWithOnePF(:,2);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOn = PlaceCellsWithOnePFOn(~isnan(PlaceCellsWithOnePFOn));
PlaceCellsWithOnePFAfter = PlaceCellsWithOnePF.PlaceCellsWithOnePF(:,3);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoAfter = PlaceCellsWithOnePFAfter(~isnan(PlaceCellsWithOnePFAfter));

sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.PlaceCellsWithOnePFinAllProt = intersect(intersect(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOff,sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOn),sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoAfter);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.PlaceCellsWithOnePFinAnyProt = union(union(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOff,sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOn),sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoAfter);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.PlaceCellsWithOnePFinOptoOfforOn = union(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOff,sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF3.OptoOn);

% saving
Folder = 'C:\MATLAB\SAVE\FinalEffect';
savePath = fullfile(Folder,sData.sessionInfo.fileID);
if ~isfolder(savePath) 
    mkdir(savePath)
end
save(savePath,'sData');
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end