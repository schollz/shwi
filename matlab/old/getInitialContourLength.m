function [startPosition] = getInitialContourLength(x,y)

    firstInd = find(x<100 & y<-50 & y>-150);
    lastInd = find(x>mean(x)+std(x) & y < 20);
    p1 = polyfit(x(firstInd),y(firstInd),1);
    try
        p2 = polyfit(x(lastInd),y(lastInd),1);
    catch
        p2(1) = 0;
        p2(2) = 0;
    end
    if abs(p2(1))>0.05
        p2(1) = 0;
        p2(2) = 0;
    end
%     p1
%     p2
    yfit1 = polyval(p1,x);
    yfit2 = polyval(p2,x);
    plot(x,y,x,yfit1,x,yfit2);
    axis([-100 500 -100 500]);
    xintersect = (p2(2)-p1(2))/(p1(1)-p2(1));
    yintersect = p1(1)*xintersect + p1(2);
    
    startPosition = xintersect;
    if abs(startPosition)>50
        startPosition = 0
    end
%     pause(1)
end