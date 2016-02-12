function [origin] = getBaseline(x,y,toPlot)
% Gets baseline, returns [xintersect, yintersect]
% [origin] = getOrigin(x,y,toPlot)

    if nargin < 3
       toPlot = 0;
    end
    addpath('../util')
    lastInd = length(y)-50:length(y);
    sortedY = sort(y);
    firstInd = find(y<sortedY(round(length(x)*0.02)) & x < mean(x)/3);
    dat = sortrows([x(firstInd) y(firstInd)]);
    p1 =  estimate_line_ver_weighted(dat(:,1),dat(:,2),1);
    
    try  
        dat = sortrows([x(lastInd) y(lastInd)]);
        p2 = polyfit(dat(:,1),dat(:,2),1);
    catch
        p2(1) = 0;
        p2(2) = 0;
    end
    if abs(p2(1))>0.05
        p2(1) = 0;
        p2(2) = 0;
    end
     yfit1 = polyval(p1,x);
     yfit2 = polyval(p2,x);
     if toPlot > 0
        plot(x,y,x,yfit1,x,yfit2,x(firstInd),y(firstInd));
        axis([-100 max(x) -300 600])
     end
    xintersect = (p2(2)-p1(2))/(p1(1)-p2(1));
    yintersect = p1(1)*xintersect + p1(2);
    origin =[xintersect,yintersect];
end