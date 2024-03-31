function newPlaceFields2(sData)
 
% Note: Grienbergers definition for new place field:
% First she looks for a Ca2+ transient that is more than 3 standard deviations (SD) about the noise (me: baseline). 
% Then, in the 6 laps after that, there needs to be a transient that is at least 3 SD about the noise in 3 out of 6 laps.
% My modified version: after induction there needs to be 5 transients at the same bin(s) in 5 out of 15 laps (30%)


mkdir(strcat(sData.sessionInfo.savePath,'\Imaging'),'NewPC2');
savePath = strcat(sData.sessionInfo.savePath,'\Imaging\NewPC2');

% set which data set you wanna use, parameters:
VelMin = 0.1;
FigVisible = 'off';
% dFF data curation: 
NoiseLevelSD = 3; % set datapoints value to zero if it is smaller, then the NoiseLevel based on SD. it was originally 3
%dFFSmoothingWindow = 1; % do a movmean smoothing on the original dFF data with bin-size specified here
%discardDFFBelow = 0.015; % discard dff below this value, seems to be unnecessary since SD level discarding is enough
RefillActivitySample = 4; % correct abrupt lack of activity, cannot handle more than 4 currently

% discard ROIs with low activity or low SNR
%RoiActivityLevel = 0.2; % discard below this
RoiSNR = 5; % discard below this
PFlengthMin = 5; % Minimal PF size in bins at induction lap
%RoiPeakdFF = 0.2; % discard below this

% set parameters for finding a newly formed place field
NoTestBinsUntil = 5; % do not test data before this bin, since the opto stimulation start around this bin. Originally = 10
TestTrialsBeforeInd = 10; % how many laps to test before induction. (should be no activity) Originally = 6
AllowedPreActivityPercentage = 4; % Out of TestTrialsBeforeInd these few trials are allowed to have activity at the same location

TestBinsBeforeIndBin = 15; % test after induction lap. how many bins to test around induction. Originally = 5
TestBinsAfterIndBin = 5; % test after induction lap. how many bins to test around induction. Originally = 5
TestTrialsAfterInd = 10; % later re-set as well, how many laps to test after induction to see activity. Originally = 20
Reliability = 70; % percentage. Originally = 50
% Shift in the place field center of mass:
TestBinsBeforeCenter = 30; % how many bins test after indution trial for the shift. Originally = 20
TestBinsAfterCenter = 10;  % how many bins test after indution trial for the shift. Originally = 10
% set
nBin = sData.behavior.meta.nBins;
nTrials = sData.behavior.wheelLapImaging-1; 
nROIs = sData.imdata.nROIs;
%signal = sData.imdata.MaoPC_deconv_IO3_RI01; % Note: I set the deconvolution parameters such, that it detects spikes if it is above 3 SNR. So every 'spike' in the deconvolution data are considered as a potential spike
signalPre = sData.imdata.roiSignals(2).dff;
signal = NaN(size(signalPre)); 
% search for transients which are 3x higher than SD of data
SD = NaN(nROIs,1);
for i = 1:1:nROIs
    %signalTemp1 = smoothdata((signalPre(i,:)),'movmean',dFFSmoothingWindow);
    signalTemp1 = signalPre(i,:);
    SD(i) = std(signalTemp1);
    signalTemp2 = signalTemp1;
    signalTemp2(signalTemp2 < NoiseLevelSD*SD(i)) = 0;
    %signalTemp2(signalTemp2 < discardDFFBelow) = 0; % discard valuse which are too small, even if larger than 2xSD. 
    signal(i,:) = signalTemp2;
end

%%% bin signal for each ROI
signalBinned = signalHeatPlot(signal,sData,VelMin,FigVisible);

% convert data to zeros (no activity) and ones (activity)
SignalZeroOne = signalBinned;
for roi = 1:1:nROIs
    SignalZeroOne{1,roi}(SignalZeroOne{1,roi} > 0) = 1; % logical transformation of the spike-matrix. 0: no activity, 1: activity
end

% Correct missing values (similar to smoothing), activity plot won't be so spotty
SignalZeroOneCorrected = SignalZeroOne;
% get rid of short lack of activity: 0 0 1 1 0 1 1 1 0 0 -> 0 0 1 1 1 1 1 1 0 0
for roi = 1:1:nROIs
    for i = 1:1:nTrials-1
        SignalZeroOneTemp = [SignalZeroOne{1,roi}(i,1:nBin) SignalZeroOne{1,roi}(i+1,1:nBin)]; %concatenate the next trial
        container = NaN(1,2);
        for j = 1:1:2*nBin-RefillActivitySample
            container(1,2) = SignalZeroOneTemp(j);
            if container(1) == 1 && container(2) == 0 && sum(SignalZeroOneTemp(j:j+RefillActivitySample))>0
                SignalZeroOneTemp(j) = 1;
                container(2) = 1;
            end
            container(1:2) = circshift(container,-1);
        end
        SignalZeroOneCorrected{1,roi}(i,1:nBin) = SignalZeroOneTemp(1:nBin); 
    end
end
%}

NewPlaceCellCollection = struct;
NewPlaceCellCollection.PF = NaN(1,10);
NewPlaceCellCollection.legend(1,1) = convertCharsToStrings('roiID'); 
NewPlaceCellCollection.legend(1,2) = convertCharsToStrings('PFindTrial'); 
NewPlaceCellCollection.legend(1,3) = convertCharsToStrings('PFindBin');
NewPlaceCellCollection.legend(1,4) = convertCharsToStrings('PFLengthBin'); 
NewPlaceCellCollection.legend(1,5) = convertCharsToStrings('PFActCenteratIndBin'); 
NewPlaceCellCollection.legend(1,6) = convertCharsToStrings('PFPeakBinatInd'); 
NewPlaceCellCollection.legend(1,7) = convertCharsToStrings('PFPeakValueatInd');
NewPlaceCellCollection.legend(1,8) = convertCharsToStrings('PFstabilaziedStartBin');
NewPlaceCellCollection.legend(1,9) = convertCharsToStrings('PFstabilaziedEndBin');
NewPlaceCellCollection.legend(1,10) = convertCharsToStrings('PFstabilaziedCenterBin');
NewPlaceCellCollection.legend(1,11) = convertCharsToStrings('PFshiftCenterBin'); 
NewPlaceCellCollection.legend(1,12) = convertCharsToStrings('veloInd');
nNewPC = 0; % number of new place fields

for roi = 1:1:nROIs
    NewPCFound = 0;
    %Controls: ROI has to be active, high enough SNR
    if sData.imdata.roiStat.signalToNoise(roi) < RoiSNR 
        continue
    end
    % generate circularized signal for each ROI
    SignalBefore = circshift(SignalZeroOneCorrected{1,roi},1,1); %SignalZeroOneSmoothed
    SignalBefore(1,:) = 0;
    SignalAfter = circshift(SignalZeroOneCorrected{1,roi},-1,1); %SignalZeroOneSmoothed
    SignalAfter(end,:) = 0;
    TempBinnedSignal = horzcat(SignalBefore,SignalZeroOneCorrected{1,roi},SignalAfter); % circularize data SignalZeroOneSmoothed
    % transform data: 0: no activity, 1: start of transient, 2: ongoing transient
    for p = 1:1:nTrials
        for q = 2:1:3*nBin
            if TempBinnedSignal(p,q)==1 && TempBinnedSignal(p,q-1)==1 || TempBinnedSignal(p,q)==1 && TempBinnedSignal(p,q-1)==2
               TempBinnedSignal(p,q) = 2;
            end
        end
    end
    
    %%% looking for new place fields 
    for i = 5:1:nTrials-TestTrialsAfterInd  % at least have 5 laps for testing pre-activity
       if NewPCFound == 1
           break % if there is new PF for the current ROI stop searching for more
       end
       for j = nBin+NoTestBinsUntil+1:1:2*nBin
          TestTrialsAfterInd = 10;
          if TempBinnedSignal(i,j) == 1 % looks for the starting point of transients, there should not be activity in the previous 6-10 laps 
               %%% Where is PF start, end, center
               PFStartTrial = i;
               PFStartBin = j-nBin;
               PFlength = 1; 
               for r = j+1:1:(4*nBin-j)
                    if TempBinnedSignal(i,r) > 0
                        PFlength = PFlength + 1;
                    else
                        break
                    end
               end
               if PFlength < PFlengthMin
                   break
               end
               PFActCenteratIndBin = floor(PFStartBin + (PFlength / 2)); % center of mass of induction PF
               if PFStartBin+PFlength > nBin
                  PFEnd = nBin;
               else
                  PFEnd =  PFStartBin+PFlength;
               end
               dFFInd = sData.imdata.binned.RoidFF{1,roi}(i,PFStartBin:PFEnd);
               PFPeakValueatInd = max(dFFInd);
               PFPeakBinatInd = find(dFFInd == max(dFFInd)) + PFStartBin -1;
               if PFStartTrial + TestTrialsAfterInd > nTrials
                   TestTrialsAfterInd = nTrials - PFStartTrial;
               end
               
               %%% test if activity on the previous 6-10 laps in the neighbouring X bins, if it was zero 
               PreActivity = zeros(TestTrialsBeforeInd,TestBinsBeforeIndBin+TestBinsAfterIndBin); 
               if i > TestTrialsBeforeInd
                  ii = i-TestTrialsBeforeInd;
               else   
                  ii = 1; % if activity was before 6-10th laps, I can test fewer laps before
               end
               for r = ii:1:i-1 % e.g. test from lap 10 to 16
                  PreActivity(r-ii+1,1:TestBinsBeforeIndBin+TestBinsAfterIndBin+1) = TempBinnedSignal(r,PFActCenteratIndBin+nBin-TestBinsAfterIndBin:PFActCenteratIndBin+nBin+TestBinsBeforeIndBin); % yes, before-after ind good as it is! 
               end
               PreActivity(PreActivity==2)=1;
               [s1,s2] = size(PreActivity);
               if sum(sum(PreActivity,2))/(s1*s2)*100 > AllowedPreActivityPercentage % if there was actity in X trials in the previous x laps, go to the next candidate
                   continue
               end
               
               % generate an array showing in which laps after induction there was activity in the same bins (+1, -1 bins) as during induction
               AfterInductionArray = TempBinnedSignal(PFStartTrial+1:(PFStartTrial+TestTrialsAfterInd),(PFActCenteratIndBin+nBin-TestBinsBeforeIndBin-1):(PFActCenteratIndBin+nBin+(TestBinsAfterIndBin-1)));
               AfterInductionArray(AfterInductionArray == 2) = 1;
               %MeanAfterIndPosTuning = mean(AfterInductionArray,1);
               AfterInductionArrayValues = sum(AfterInductionArray,2); % can be estimated the amplitude of activation if needed
               AfterInductionArraySum = max(AfterInductionArray,[],2); 
               if sum(AfterInductionArraySum) >= TestTrialsAfterInd*(Reliability/100) && sum(AfterInductionArrayValues) >= (TestBinsBeforeIndBin + TestBinsAfterIndBin)*TestTrialsAfterInd/5
                  % considered as newly formed place field, if x% of the induction following trials have transients, plus the activation is large enough
                  nNewPC = nNewPC + 1;
                  NewPlaceCellCollection.PF(nNewPC,1) = roi; % write the ROI number into column one is it is a place cell
                  NewPlaceCellCollection.PF(nNewPC,2) = PFStartTrial; % induction lap
                  NewPlaceCellCollection.PF(nNewPC,3) = PFStartBin; % induction bin
                  NewPlaceCellCollection.PF(nNewPC,4) = PFlength; % induction transient length (bin) 
                  NewPlaceCellCollection.PF(nNewPC,5) = PFActCenteratIndBin; % Center of Mass of induction lap transient
                  NewPlaceCellCollection.PF(nNewPC,6) = PFPeakBinatInd;
                  NewPlaceCellCollection.PF(nNewPC,7) = PFPeakValueatInd;
                  % calculate center of mass and length of established place field after the induction laps
                  %PFstabilaziedStartBinArray = NaN(TestTrialsAfterInd,1);
                  %PFstabilaziedEndBinArray = NaN(TestTrialsAfterInd,1);
                  % calculate the center of mass of activity
                  TempAct = TempBinnedSignal(PFStartTrial+1:PFStartTrial+TestTrialsAfterInd,nBin+PFActCenteratIndBin-TestBinsBeforeCenter:nBin+PFActCenteratIndBin+TestBinsAfterCenter);                      
                  TempAct(TempAct == 2) = 1;
                  PosTuningForShift = mean(TempAct,1);
                  PFstabilizedCenterBin = (PFActCenteratIndBin-TestBinsBeforeCenter-1) + floor(mean(find(PosTuningForShift==max(PosTuningForShift))));
                  PFshiftCenterBin =  PFstabilizedCenterBin - PFActCenteratIndBin;
                  %{
                  % find starting bin:
                          for h = 2:1:TestBinsBeforeCenter+TestBinsAfterCenter
                                if ActivityTemp(h) == 1 && ActivityTemp(h-1) == 0 % Where is stabilized PF start
                                    PFstabilaziedStartBin = PFActCenteratIndBin - TestBinsBeforeCenter + h - 1;
                                end
                                if ActivityTemp(h) > 0 && ActivityTemp(h+1) == 0 % Where is stabilized PF end
                                    PFstabilaziedEndBin = PFActCenteratIndBin - TestBinsBeforeCenter + h;
                                    PFstabilaziedStartBinArray(n) = PFstabilaziedStartBin;
                                    PFstabilaziedEndBinArray(n) = PFstabilaziedEndBin;
                                    %PFshiftCenterBinArray(n) = round(PFstabilaziedStartBin+((PFstabilaziedEndBin - PFstabilaziedStartBin)/2));
                                    %PFestabLengthArray(n) = PFstabilaziedEndBin - PFstabilaziedStartBin;
                                end                        
                          end
                  %}                      
                  
                  %NewPlaceCellCollection.PF(nNewPC,8) = nanmean(PFstabilaziedStartBinArray);
                  %NewPlaceCellCollection.PF(nNewPC,9) = nanmean(PFstabilaziedEndBinArray);
                  NewPlaceCellCollection.PF(nNewPC,10) = PFstabilizedCenterBin;
                  NewPlaceCellCollection.PF(nNewPC,11) = PFshiftCenterBin;
                  NewPlaceCellCollection.PF(nNewPC,12) = sData.behavior.binning.veloBinned(i,PFPeakBinatInd);
                  NewPCFound = 1;
                  break
               end
          end
       end
    end
end

% check if induction lap was in opto-off, opto-on, or opto-after trial
OptoOn = 0; OptoOff = 0; OptoAfter = 0;
for i = 1:1:nNewPC
    if ismember(NewPlaceCellCollection.PF(i,2),sData.behavior.opto.OptoOnTrialsIndices) == 1
        OptoOn = OptoOn + 1;
    elseif ismember(NewPlaceCellCollection.PF(i,2),sData.behavior.opto.OptoOffTrialsIndices) == 1
        OptoOff = OptoOff + 1;
    elseif ismember(NewPlaceCellCollection.PF(i,2),sData.behavior.opto.AfterOptoTrialsIndices) == 1
        OptoAfter = OptoAfter + 1;
    end
    NewPlaceCellCollection.InductionTrial.OptoOn = OptoOn;
    NewPlaceCellCollection.InductionTrial.OptoOff = OptoOff;   
    NewPlaceCellCollection.InductionTrial.OptoAfter = OptoAfter;
    NewPlaceCellCollection.InductionTrial.OptoOnOffRatio = OptoOn/OptoOff;
    NewPlaceCellCollection.InductionTrial.OptoOnOffRatioNormalized = (OptoOn/length(sData.behavior.opto.OptoOnTrialsIndices)) / (OptoOff/length(sData.behavior.opto.OptoOffTrialsIndices));
end

%%% figures
%{
%%% Plot current ROI activity
for i=1:1:nNewPC
    roi = NewPlaceCellCollection.PF(i,1)
    DataForPlotting = SignalZeroOneCorrected{1, roi};
    %DataForPlotting(DataForPlotting==2) = 1;
    figure('Color','white')
    imagesc(DataForPlotting(:,1:nBin))
    title(strcat('ROI-',num2str(roi)));
    FileName = strcat(sData.sessionInfo.fileID,'NewPCForm-ROI',num2str(roi));
    %savefig(fullfile(savePath,FileName));
    saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));
end
    close all; 

% induction velo vs PF width 
figure('Color','white');
scatter(NewPlaceCellCollection.PF(:,7),NewPlaceCellCollection.PF(:,8)*sData.behavior.meta.binSize)
xlabel('Velocity during induction (cm/s)');
ylabel('Place Field Width (cm)');
title(strcat(sData.sessionInfo.fileID,'-Velocity vs established place field width'));
FileName = strcat(sData.sessionInfo.fileID,'-VeloVsPFwidth');
%savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% PFwidth shift 
figure('Color','white');
histogram(NewPlaceCellCollection.PF(:,6)*sData.behavior.meta.binSize)
xlabel('Shift in place tuning after induction (cm)');
ylabel('Count (place fields)');
title(strcat(sData.sessionInfo.fileID,'-Shift In Place Tuning'));
FileName = strcat(sData.sessionInfo.fileID,'-ShiftInPlaceTuning');
%savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Cumulative distribution plot for place field formation
% convert nTrials into a percentage where 100% is the end of session
MaxTrials = nTrials - TestTrialsAfterInd;
CumDistr = NewPlaceCellCollection.PF(:,2)/MaxTrials*100;

figure('Color','white');
cdfplot(CumDistr) 
xlabel('Lap of place field induction');
ylabel('Cumulative fraction of place cells');
title(strcat(sData.sessionInfo.fileID,'-Cumulative distribution of lap of place field induction'));
FileName = strcat(sData.sessionInfo.fileID,'-CumDistPF');
%savefig(fullfile(savePath,FileName));
saveas(gcf,(fullfile(savePath,[FileName '.jpg'])));

% Save file to same path where LV files can be found 
save(fullfile(savePath,strcat(sData.sessionInfo.fileID,'_NewPlaceCellCollection.mat')),'NewPlaceCellCollection');
%}

end


