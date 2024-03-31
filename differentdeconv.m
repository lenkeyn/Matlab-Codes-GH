
% create an array where the first column is the ROI-ID, the second is the SNR with which the deconvolution should be made    
[nROIs, nSamples] = size(dff);
RoiSnrArrayPre = NaN(nROIs,2); 
% load data from excel table wher SNR values are set !!!! do manually!!!!
RoiSnrArray = RoiSnrArrayPre(~isnan(RoiSnrArrayPre(:,2)),1:2);

cia_dec = zeros(size(dff));
cia_den = zeros(size(dff));
cia_opt{1,nROIs} = struct;
for i = 1:1:length(RoiSnrArray)
   roi = RoiSnrArray(i,1);
   snr = RoiSnrArray(i,2);
   opt = struct('lam_pr', 0.5, 'spk_SNR', snr, 'fps', 31, 'type', 'ar1', 'tau_dr', [550, 180]./ 1000 .* opt.fps);
   dff_temp = dff(roi,:);
   [cia_dec_temp, cia_den_temp, cia_opt_temp] = signalExtraction.CaImAn.deconvolve(dff_temp, opt);
   %cia_dec(roi, :) = cia_dec_temp;
   %cia_den(roi, :) = cia_den_temp;
   %cia_opt{1,roi} = cia_opt_temp;
end
 % check a roi dff and deconvolved data
 roi = 2;
 figure();
 plot(dff(roi,:)); hold on
 plot(cia_dec(roi,:))
 
cia_dec(roi,10001:25120) = cia_dec_temp(10001:25120);
cia_den(roi,10001:25120) = cia_den_temp(10001:25120);
cia_opt{2,roi} = cia_opt_temp;
 
 
 
 
