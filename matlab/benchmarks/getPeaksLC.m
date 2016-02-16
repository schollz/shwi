function [pks,locs] = getPeaksSEGM(xbaseline,Fbaseline,toPlot)
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
xbaseline=data(:,1)
Fbaseline=data(:,2)
toPlot = 1;

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
[PKS,LOCS] = findpeaks(f,xi,'MinPeakProminence',0.0001,'MinPeakHeight',1e-3,'MinPeakDistance',10);
L = unique(LOCS);
F = PKS;

Fpks = [];
Xpks = [];
for i=1:length(L)
    yMax = WLC4(0.3,L(i),x+3)+3;
    yMin = WLC4(0.5,L(i),x-3)-3;
%     plot(x,y,x,yMin,x,yMax)
%     axis([0 max(x) -50 max(y)])
%     pause(1)
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
    plot(xbaseline,Fbaseline,locs,pks,'or')
    subplot(2,1,2)
    plot(xi,f,LOCS,PKS,'o')
end


end