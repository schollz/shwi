function [pks,locs,allLcs] = getPeaksLC(xbaseline,Fbaseline,toPlot)
% Gets peaks using the SEGM method
% [pks,locs] = getPeaksSEGM(xbaseline,Fbaseline,toPlot)
    if nargin < 3
       toPlot = 0;
    end
%     
% 
% xbaseline = r{1}.x;
% Fbaseline= r{1}.y;
% toPlot = 1;
% xbaseline=data(:,1)
% Fbaseline=data(:,2)
% toPlot = 1;

x = xbaseline;
y = Fbaseline;

allLcs = [];
for di=1:length(y)
    if y(di) > 10
        try
            LC0 = getLc2(0.4,x(di),y(di));
            allLcs = [allLcs; LC0];
        catch
        end
    end
end
[f,xi]=ksdensity(allLcs,linspace(0,max(x)*1.5,length(x)),'bandwidth',4);
[PKS,LOCS] = findpeaks(f,xi,'MinPeakProminence',0.0005,'MinPeakHeight',4e-3,'MinPeakDistance',10);
L = unique(LOCS);
F = PKS;

Fpks = [];
Xpks = [];
for i=1:length(L)
    yMax = WLC4(0.3,L(i),x+3)+3;
    yMin = WLC4(0.5,L(i),x-3)-3;

%     subplot(2,1,1)
%     hold on
%     plot(x,y,'k',x(find(x<L(i)*0.9)),yMin(find(x<L(i)*0.9)),'r',x(find(x<L(i)*0.9)) ,yMax(find(x<L(i)*0.9)),'r')
%     axis([-10 max(x)*1.05 -50 max(y)*1.1])
%     pause(1)
%     title('Force-extension')
%     legend('Data','Peaks','Lc-space limits')
%     xlabel('Extension [nm]')
%     ylabel('Force [pN]')
%     subplot(2,1,2)
%     axis([-10 max(x)*1.05 0 0.01])
%     ylabel('P')
%     title('Contour-length smoothed histogram')
%     xlabel('Lc [nm]')
    
    yInd = find(y>yMin & y < yMax);
    [newY,ind] = max(y(yInd));
    newX = x(yInd(ind));
    if sum(Xpks==newX) == 0
        if i > 1
            if abs(newX-Xpks(end)) < 10
            else
                Xpks = [Xpks; newX];
                Fpks = [Fpks; newY];            
            end
        else
            Xpks = [Xpks; newX];
            Fpks = [Fpks; newY];
        end
    end
end


locs = Xpks;
pks = Fpks;
if toPlot > 0
    subplot(2,1,1)
    plot(xbaseline,Fbaseline,'k',locs,pks,'or')
    axis([-10 max(x)*1.05 -50 max(y)*1.1])
    title('Force-extension')
    xlabel('Extension [nm]')
    ylabel('Force [pN]')
    subplot(2,1,2)
    plot(xi,f,LOCS,PKS,'o')
    axis([-10 max(x)*1.05 0 max(f)*1.1])
    ylabel('P')
    title('Contour-length smoothed histogram')
    xlabel('Lc [nm]')
end


end