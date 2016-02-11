function [pks,locs] = findFEpeaks(xbaseline,Fbaseline)
    xbaseline = xbaseline(find(Fbaseline>-50));
    Fbaseline = Fbaseline(find(Fbaseline>-50));


    y=Fbaseline;
    u = ones(size(y));
    bestSegm = 0;
    bestNum = 100000;
    for R2 = 5:3:20
    segm = segment([y,u],[0 1 1],R2);
    if sum(isnan(segm))<bestNum
        bestNum = sum(isnan(segm));
        bestSegm = segm;
    end
    if sum(isnan(segm)) < 500
        break
    end
    end
    segm = bestSegm;

    currentNum = segm(1);
    currentXs = [];
    newX = [];
    newY = [];
    for i=1:length(segm)
       if segm(i)==currentNum && isnan(segm(i))==0
           currentXs = [currentXs; xbaseline(i)];
       else
           if length(currentXs) > 1
               newX = [newX; mean(currentXs)];
               newY = [newY; currentNum];
           end
           currentXs = [xbaseline(i)];
           currentNum = segm(i);
       end
    end
    newX(isnan(newX)) = 0;
    dats = [newX+rand(size(newX))/100 newY];
    dats = sortrows(dats);
    newX = dats(:,1);
    newY = dats(:,2);

    [pks,locs] = findpeaks(newY,newX,'MinPeakProminence',6,'MinPeakDistance',12,'MinPeakHeight',11);
%     subplot(2,1,1)
%     plot(xbaseline,Fbaseline,xbaseline,segm,'.',locs,pks,'or',newX,newY,'m*-')
%     subplot(2,1,2)
%     plot(locs,pks,'or',newX,newY,'m-')
    
%      axis([-100 500 -150 250])


end