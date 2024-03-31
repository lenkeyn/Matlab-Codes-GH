function sData = LandmarkCellDetectionOptoOnDeconvVIPpaper(sData)

%%% CHECK FOR LANDMARK CELLS: 
% there are two set of landmarks on the wheel at certain distance: 
% velcro at 26 and 86 cm and hot glue spikes at 46 and 66 cm. All tactile cues span for 2 cms.
% The distance between HGS are 20 cm and velcros are 60 cm. I look for double peaks with these distances. 
% Since in the dF/F data the second peak is often masked by the first peak's decay, I check these double peaks in the dconvolved (unsmoothed) data

% I will check landmark properties in each opto protocol (off, on, after)

LandmarkArea = floor(10/sData.behavior.meta.binSize); % 10 cm (changed from 4 cm), search for landmark cells where the landmark was presented and a little addition
VelcroStart1 = floor(26/sData.behavior.meta.binSize); % rounding using floor instead of ceil , for all other
VelcroEnd1 = VelcroStart1 + LandmarkArea;
VelcroStart2 = floor(86/sData.behavior.meta.binSize);
VelcroEnd2 = VelcroStart2 + LandmarkArea;
HGSStart1 = floor(46/sData.behavior.meta.binSize);
HGSEnd1 = HGSStart1 + LandmarkArea;
HGSStart2 = ceil(66/sData.behavior.meta.binSize);
HGSEnd2 = HGSStart2 + LandmarkArea +1; % delay in Ca transient due to previous Ca transient evoked by first HGS
DiffVelcro = floor(60/sData.behavior.meta.binSize); %30 bin
DiffHGS = floor(20/sData.behavior.meta.binSize); %10 bin
AllowedJitter = 1; % bins
MinPeakDetectionMainPeak = 0.001; % spike detection treshold in the deconvolved dataset, it was 0.0001
MinPeakDetectionSecondPeak = 0.0003; % spike detection treshold in the deconvolved dataset, it was 0.0001

datasetDeconv = sData.imdata.binned.RoiDeconvolved;
LandmarkDetection2 = struct;

LandmarkCells2 = struct;    
LandmarkCells2.LandmarkCells = NaN(sData.imdata.nROIs,3); 
LandmarkCells2.VelcroCells = NaN(sData.imdata.nROIs,3);
LandmarkCells2.HGSCells = NaN(sData.imdata.nROIs,3);
LandmarkCells2.Peak1PosBin = NaN(sData.imdata.nROIs,3);
LandmarkCells2.Peak2PosBin = NaN(sData.imdata.nROIs,3);

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
            DataTemp1 = [nanmean(datasetDeconv{1,i}(Trials,:),1) nanmean(datasetDeconv{1,i}(Trials,:),1)];  % concatentate the data twice to be able to use as circular dataset
            [pks,locs] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionMainPeak); % peaks and indices at which the peaks occur. Very local peaks are also detected
            PeakAmplitudeLocation = locs(find(pks == max(pks),1)); % search for the peak amplitude's location (which bin)
            if isempty(PeakAmplitudeLocation)
                continue
            end
            % check if it can be a velcro or a HGS cell, then check if there is a peak at the other similar landmark (one bin jitter allowed for Velcro, 2 bin for HGS)
            if PeakAmplitudeLocation >= VelcroStart1 &&  PeakAmplitudeLocation <= VelcroEnd1 % Velcro1
                DiffVelcroArray = PeakAmplitudeLocation + DiffVelcro - AllowedJitter : PeakAmplitudeLocation + DiffVelcro + AllowedJitter;
                DiffVelcroArray(DiffVelcroArray < VelcroStart2) = NaN;
                DiffVelcroArray(DiffVelcroArray > VelcroEnd2) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if  any(ismember(locs2,DiffVelcroArray)>0)
                    LandmarkCells2.LandmarkCells(i,k) = i; 
                    LandmarkCells2.VelcroCells(i,k) = i;
                    LandmarkCells2.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells2.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffVelcroArray)>0));
                end
            elseif PeakAmplitudeLocation >= VelcroStart2 &&  PeakAmplitudeLocation <= VelcroEnd2 %Velcro2
                DiffVelcroArray = PeakAmplitudeLocation - DiffVelcro - AllowedJitter : PeakAmplitudeLocation - DiffVelcro + AllowedJitter;
                DiffVelcroArray(DiffVelcroArray < VelcroStart1) = NaN;
                DiffVelcroArray(DiffVelcroArray > VelcroEnd1) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffVelcroArray)>0) 
                    LandmarkCells2.LandmarkCells(i,k) = i; 
                    LandmarkCells2.VelcroCells(i,k) = i;
                    LandmarkCells2.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells2.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffVelcroArray)>0));
                end
            elseif PeakAmplitudeLocation >= HGSStart1 &&  PeakAmplitudeLocation <= HGSEnd1
                DiffHGSArray = PeakAmplitudeLocation + DiffHGS - AllowedJitter : PeakAmplitudeLocation + DiffHGS + AllowedJitter + 1; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
                DiffHGSArray(DiffHGSArray < HGSStart2) = NaN;
                DiffHGSArray(DiffHGSArray > HGSEnd2) = NaN;
                [Peaks2,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffHGSArray)>0) 
                    LandmarkCells2.LandmarkCells(i,k) = i; 
                    LandmarkCells2.HGSCells(i,k) = i;
                    LandmarkCells2.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells2.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffHGSArray)>0));
                end
            elseif PeakAmplitudeLocation >= HGSStart2 &&  PeakAmplitudeLocation <= HGSEnd2
                DiffHGSArray = PeakAmplitudeLocation - DiffHGS - AllowedJitter -1 : PeakAmplitudeLocation - DiffHGS + AllowedJitter; % an extra bin is added for jitter because of the Ca transient of the previous HGS decays here and mask new transient
                DiffHGSArray(DiffHGSArray < HGSStart1) = NaN;
                DiffHGSArray(DiffHGSArray > HGSEnd1) = NaN;
                [~,locs2] = findpeaks(DataTemp1,'MinPeakProminence',MinPeakDetectionSecondPeak);
                if any(ismember(locs2,DiffHGSArray)>0)
                    LandmarkCells2.LandmarkCells(i,k) = i; 
                    LandmarkCells2.HGSCells(i,k) = i;
                    LandmarkCells2.Peak1PosBin(i,k) = PeakAmplitudeLocation;
                    LandmarkCells2.Peak2PosBin(i,k) = min(locs2(ismember(locs2,DiffHGSArray)>0));
                end
            end
        end
    end     
end

for m = 1:1:3
    LandmarkCellListPre = LandmarkCells2.LandmarkCells(1:end,m);
    VelcroCellsListPre =  LandmarkCells2.VelcroCells(1:end,m);
    HGSCellsListPre = LandmarkCells2.HGSCells(1:end,m);
    Peak1PosBinPre = LandmarkCells2.Peak1PosBin(1:end,m);
    Peak2PosBinPre = LandmarkCells2.Peak2PosBin(1:end,m);
    if m == 1
        LandmarkDetection2.OptoOff.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetection2.OptoOff.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetection2.OptoOff.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetection2.OptoOff.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetection2.OptoOff.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    elseif m == 2
        LandmarkDetection2.OptoOn.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetection2.OptoOn.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetection2.OptoOn.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetection2.OptoOn.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetection2.OptoOn.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    elseif m == 3
        LandmarkDetection2.OptoAfter.LandmarkCellsList = LandmarkCellListPre(~isnan(LandmarkCellListPre)); 
        LandmarkDetection2.OptoAfter.Peak1PosBin = Peak1PosBinPre(~isnan(Peak1PosBinPre));
        LandmarkDetection2.OptoAfter.Peak2PosBin = Peak2PosBinPre(~isnan(Peak2PosBinPre));
        LandmarkDetection2.OptoAfter.VelcroCellsList = VelcroCellsListPre(~isnan(VelcroCellsListPre)); 
        LandmarkDetection2.OptoAfter.HGSCellsList = HGSCellsListPre(~isnan(HGSCellsListPre));
    end
end
       

LandmarkDetection2.Params.LandmarkArea = LandmarkArea; 
LandmarkDetection2.Params.VelcroStart1 = VelcroStart1;
LandmarkDetection2.Params.VelcroEnd1 = VelcroEnd1;
LandmarkDetection2.Params.VelcroStart2 = VelcroStart2;
LandmarkDetection2.Params.VelcroEnd2 = VelcroEnd2;
LandmarkDetection2.Params.HGSStart1 = HGSStart1;
LandmarkDetection2.Params.HGSEnd1 = HGSEnd1;
LandmarkDetection2.Params.HGSStart2 = HGSStart2;
LandmarkDetection2.Params.HGSEnd2 = HGSEnd2;
LandmarkDetection2.Params.DiffVelcro = DiffVelcro;
LandmarkDetection2.Params.DiffHGS = DiffHGS;
LandmarkDetection2.Params.AllowedJitter = AllowedJitter; 
LandmarkDetection2.Params.MinPeakDetection = MinPeakDetectionSecondPeak;


LandmarkDetection2.LandmarkCellinAllProt = intersect(intersect(LandmarkDetection2.OptoOff.LandmarkCellsList,LandmarkDetection2.OptoOn.LandmarkCellsList),LandmarkDetection2.OptoAfter.LandmarkCellsList);%,MaoOpto.OptoAfter.PlaceCells);
LandmarkDetection2.LandmarkCellinAnyProt = union(union(LandmarkDetection2.OptoOff.LandmarkCellsList,LandmarkDetection2.OptoOn.LandmarkCellsList),LandmarkDetection2.OptoAfter.LandmarkCellsList);%,MaoOpto.OptoAfter.PlaceCells);
LandmarkDetection2.LandmarkCellinOptoOfforOptoOn = union(LandmarkDetection2.OptoOff.LandmarkCellsList,LandmarkDetection2.OptoOn.LandmarkCellsList);
sData.imdata.MaoPC_Opto_dff.LandmarkCells2 = LandmarkDetection2;
        

%%% detect place cells with only one place field
datasetDFF = sData.imdata.binned.RoidFF; % use dFF data

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
    
    PlaceCellsWithOnePF = struct;
    PlaceCellsWithOnePF.PlaceCellsWithOnePF = NaN(sData.imdata.nROIs,1);
    PlaceCellsWithOnePF.PeakPosdFFBin = NaN(sData.imdata.nROIs,1);
    PlaceCellsWithOnePF.PeakPosDeconvBin = NaN(sData.imdata.nROIs,1);
    
    for i = 1:1:sData.imdata.nROIs % search for place cells with one peak 
        if PCdff.Criteria123Passed(i) > 0 && ~any(LandmarkDetection2.LandmarkCellinAnyProt==i)  % if the cell was a place cell but not a landmark cell
            %DataTemp = smoothdata(nanmean(dataset{1,i},1),'Gaussian',3);  % use little bit of smoothing on deconvolved data, because sometimes larger transients cause double peaks in deconvolution
            %[pks,locs] = findpeaks(DataTemp,'MinPeakProminence',MinPeakDetection); % peaks and indices at which the peaks occur.
            DataTemp2 = [nanmean(datasetDFF{1,i}(Trials,:),1) nanmean(datasetDFF{1,i}(Trials,:),1)];  % use little bit of smoothing on deconvolved data, because sometimes larger transients cause double peaks in deconvolution
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
                PlaceCellsWithOnePF.PlaceCellsWithOnePF(i) = i;
                PlaceCellsWithOnePF.PeakPosdFFBin(i) = PeakAmplitudeLocation;
                if PeakAmplitudeLocation > sData.behavior.meta.nBins
                    PlaceCellsWithOnePF.PeakPosdFFBin(i) = PeakAmplitudeLocation - sData.behavior.meta.nBins;
                end
                % find the peak location in the deconvolved signal
                DataTemp3 = nanmean(datasetDeconv{1,i},1);
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
                PlaceCellsWithOnePF.PeakPosDeconvBin(i) = locs2(locs3 == min(locs3));
            end
        end
    end
    PlaceCellsWithOnePF.PlaceCellsWithOnePFList = PlaceCellsWithOnePF.PlaceCellsWithOnePF(~isnan(PlaceCellsWithOnePF.PlaceCellsWithOnePF));

    % write data into struct
    for m = 1:1:3
        if m == 1
            sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff = PlaceCellsWithOnePF;
        elseif m == 2
            sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOn = PlaceCellsWithOnePF;
        elseif m == 3
            sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoAfter = PlaceCellsWithOnePF;
        end
    end
    
end

sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.PlaceCellsWithOnePFinAllProt = intersect(intersect(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff.PlaceCellsWithOnePFList,sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOn.PlaceCellsWithOnePFList),sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoAfter.PlaceCellsWithOnePFList);
sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.PlaceCellsWithOnePFinAnyProt = union(union(sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOff.PlaceCellsWithOnePFList,sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoOn.PlaceCellsWithOnePFList),sData.imdata.MaoPC_Opto_dff.PlaceCellsWithOnePF.OptoAfter.PlaceCellsWithOnePFList);


% saving
% save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');

end