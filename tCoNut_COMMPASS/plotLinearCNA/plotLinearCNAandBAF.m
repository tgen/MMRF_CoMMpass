function plotLinearCNAandBAF(cnaTSV,bafTXT,outName)
%plotLinearCNAandBAF(cnaTSV,bafTXT,outName)
% Defaults for this blog post
width = 10;     % Width in inches
height = 5;    % Height in inches
alw = 0.75;    % AxesLineWidth
fsz = 11;      % Fontsize
lw = 1.5;      % LineWidth
msz = 8;       % MarkerSize

% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*100]);

% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close(gcf);

cgh = readTSV(cnaTSV);
baf = readBAF(bafTXT);
bafFreq = baf(:,[1 2 11]);

endChr=max(bafFreq(:,1));

hFig = figure;
set (hFig,'color','w');
subplot(2,1,1)
plotLinearBAF(bafFreq(:,3),bafFreq(:,1),unique(bafFreq(:,1))','B-Allele Frequencies (BAF)',[0 0.4470 0.7410])
subplot(2,1,2)
plotLinearCGH(cgh(:,3),cgh(:,1),unique(cgh(:,1))','Copy Number',[0.4660 0.6740 0.1880])


print(hFig,outName,'-dpng','-r300');
close(hFig);
end










