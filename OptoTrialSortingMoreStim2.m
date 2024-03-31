function sData = OptoTrialSortingMoreStim2(sData,DiscardCmBeginning)

if any(ismember(fields(sData.daqdata),'optoStimProtocol'))

    % downsample optoStimProtocol data, which shows which protocol (trial) types were used during the experiment. 
    % Usually protocol zero is the control (later I label control as protocol 1), protocol 1 is the highest intensity (if there are more). For details of protocols see sData.stimProtocols
    optoStimDS = sData.daqdata.optoStimProtocol(sData.behavior.details.frameStartIndices);
    if ~any(optoStimDS==1) % sometimes I wrongly set protocol 1 as protocol 2 during recording, when there is only one protocol, correction
        optoStimDS(optoStimDS==2)=1;
        optoStimDS(optoStimDS==3)=2;
        optoStimDS(optoStimDS==4)=3;
    end
    sData.behavior.optoStimProtocolDS = optoStimDS;
    nTrials = sData.behavior.wheelLapImaging-1;
    nBins = sData.behavior.meta.nBins;
    %DiscardCmBeginning = 10; % discard first 10 cm, because many times full trial stimulation shifts to beginning of next lap

    % binnning of optoStim data
    OptoStimProtMatrix = NaN(nTrials,nBins); % OptoStimProtMatrix : what protocol was set. Real optical stimulation: sData.behavior.opto.OptoStimOnMatrix
    OptoStimProtMatrixReal = NaN(nTrials,nBins);
    % calculate the applied stim protocol during each bin 
    for i = 1:1:nTrials  % rows are trials
        for j = 1:1:nBins  % columns (distance bins)  
            OptoStimProtMatrix(i,j) = mean(optoStimDS(sData.behavior.binning.enterIntoBinIndexExtended(i,j):(sData.behavior.binning.enterIntoBinIndexExtended(i,j)+round(sData.behavior.binning.SampleSpentInBinExtendedBins(i,j))-1))); % what was set as a protocol
            OptoStimProtMatrixReal(i,j) = mean(optoStimDS(sData.behavior.binning.enterIntoBinIndexExtended(i,j):(sData.behavior.binning.enterIntoBinIndexExtended(i,j)+round(sData.behavior.binning.SampleSpentInBinExtendedBins(i,j))-1)));
            if OptoStimProtMatrixReal(i,j)>0 && sData.behavior.opto.OptoStimOnMatrix(i,j) == 0 % if there were no stimulation indeed , despite the protocol was set to be
                OptoStimProtMatrixReal(i,j) = 0;
            end
        end
    end
    for i = 1:1:nTrials
        if sData.behavior.opto.FailedOptoTrials(i) == 1 % if the stimulation failed , discard the whole trial
            OptoStimProtMatrixReal(i,:) = -1; %!!! changed from NaN to 0 at 2021.06.05, -1 on 2021.08.05. , means that that trial was meant to be optically stimulated, but failed
            %OptoStimProtMatrixReal(i,OptoStimProtMatrixReal(i,:)>0) = -1; % code was until 2023.01.19.
        end
        if sData.behavior.opto.AfterOptoTrials(i) == 1 % separate after-opto trials
            OptoStimProtMatrixReal(i,:) = 4; %!!! changed from NaN to 0 at 2021.06.05, -1 on 2021.08.05. , means that that trial was meant to be optically stimulated, but failed
        end
    end
    StartingBin = round(DiscardCmBeginning/sData.behavior.meta.binSize)+1; % do not count the first e.g.10 cms (in bins), because many time stimulation from previous trial ends
    OptoStimProtTrialsWhatWasSet = round(nanmean(OptoStimProtMatrix(1:nTrials,StartingBin:nBins),2))+1; % discard the first X bins, since sometimes stimulation remained from previous trial.  I added +1 to start the protocol number with 1, not zero
    OptoStimProtTrialsReal = round(nanmean(OptoStimProtMatrixReal(1:nTrials,StartingBin:nBins),2))+1; % discard the X three bins, since sometimes stimulation remained from previous trial.  I added +1 to start the protocol number with 1, not zero
    % -1: after opto, 0: failed trial, 1: ctr, 2-3-4: optical stimulation
    
    sData.behavior.optoMoreProts.OptoStimProtMatrixWhatWasSet = OptoStimProtMatrix;
    sData.behavior.optoMoreProts.OptoStimProtTrialsWhatWasSet = OptoStimProtTrialsWhatWasSet;
    
    sData.behavior.optoMoreProts.OptoStimProtMatrixReal =  OptoStimProtMatrixReal;
    sData.behavior.optoMoreProts.OptoStimProtTrialsReal = OptoStimProtTrialsReal;
    


    % PLOT FIGURE
    % theoretical optical stimulation, which was set 
    figure('Color','white'); 
    imagesc(1:160,1:nTrials,OptoStimProtMatrix) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    mymap = [
         0.5 0.5 0.5 % failed trials dark grey
         0.8 0.8 0.8 % control grey
         0.9 0.2 0.2 % opto protocol 1 red
         1 0.5 0 % opto protocol 2 orange
         0.3 0.3 0.3
         1 1 0 % opto protocol 3 yellow
         ]; % set colormap colors, 
    colormap(mymap);
    xlabel('Position on wheel (cm)');
    ax = gca; ax.TickDir = 'out';
    xticks([0,25,50,75,100,125,150]);
    ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID  ,'-optical-protocols-WhatWasSet'));
    FileName = strcat('OptoStimProtWhatWasSet-',sData.sessionInfo.fileID);
    savefig(fullfile(sData.sessionInfo.savePath,FileName));
    saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));
    
    % real optical stimulation, which was actually happened 
    figure('Color','white'); 
    imagesc(1:160,1:nTrials,OptoStimProtMatrixReal) %(1:number of bins;1:number of trials)   imagesc(1:CR,1:(TR-1),(BinLick))
    mymap = [
         0.6 0.6 0.6 % 1 control grey
         0.8 0.8 0.8 % 5 after trials light gray
         0.9 0.2 0.2 % 2 opto protocol 1 red
         1 1 0 % 3 opto protocol 2 yellow
         1 0.5 0 % 4 opto protocol 3 orange
         0 0 0 %fail black?
         ]; % set colormap colors, 
    colormap(mymap);
    xlabel('Position on wheel (cm)');
    ax = gca; ax.TickDir = 'out';
    xticks([0,25,50,75,100,125,150]);
    ylabel('Trials');
    title(strcat(sData.sessionInfo.fileID  ,'-optical-protocols-real'));
    FileName = strcat('OptoStimProtReal-',sData.sessionInfo.fileID);
    savefig(fullfile(sData.sessionInfo.savePath,FileName));
    saveas(gcf,(fullfile(sData.sessionInfo.savePath,[FileName '.jpg'])));
    
    % Save file to same path where other files can be found 
    save(fullfile(sData.sessionInfo.savePath,strcat(sData.sessionInfo.fileID,'_sData.mat')),'sData');
    
end

end
