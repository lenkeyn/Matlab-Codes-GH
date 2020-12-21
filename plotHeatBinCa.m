function plotHeatBinCa(data,fileID,roi,ylab,BinSize,BinNu,TRNu,visibility) 

% parameters: x=160; y=CADATA.TRNu; data=CADATA.VeloLimInBin;
% data: matrix in which rows are trials, columns are bins, binned velocity values in cells
% ylabel : if Velo set it to 'Speed (cm/s)' , if dff set it to 'dF/F'
% Xaxis = 1:2:155,Yaxis = 1:TRNu

%parts = strsplit(filePath, '\');
%filedir = parts{end-2};
if roi == 0
    figtitle = fileID;
else
    figtitle = strcat(fileID,sprintf(' ROI #%d',roi));
end


Xaxis = 0:BinSize:BinSize*BinNu;
%PLOT FIGURE
figure('Color','white','visible',visibility); 
%plot:
%imagesc(1:160,1:CADATA.TRNu,(CADATA.VeloLimInBin)) 
%imagesc(1:x,1:y,data) %(1:number of bins;1:number of trials)
imagesc(Xaxis,1:TRNu,data);
c = colorbar;
colormap(jet);
c.Label.String = ylab;
c.Label.FontSize = 11;
c.TickDirection = 'out'; 
%if ylab == 'Licks/cm'
%   caxis([0 5]);
%end    
%caxis([CADATA.VelMin 50]); %set limits for color plot, below 1st black, above 2nd white
axis([0 157 0 TRNu]);
xlabel('Position on Wheel (cm)');
ax = gca;
ax.TickDir = 'out';
xticks([0,25,50,75,100,125,150]);
ylabel('Trials');
title(figtitle);

end

%{
mymap = [0 0 0; 0.7 0 0; 0.8 0 0; 0.9 0 0; 1 0 0;
    1 0.1 0; 1 0.2 0; 1 0.3 0; 1 0.4 0; 1 0.5 0; 1 0.6 0; 1 0.7 0; 1 0.8 0; 1 0.9 0; 1 1 0;
    1 1 0.2; 1 1 0.4; 1 1 0.6; 1 1 0.8];
%}