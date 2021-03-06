function [pks,locs] = getPeaksSEGM(xbaseline,Fbaseline,toPlot)
% Gets peaks using the SEGM method
% [pks,locs] = getPeaksSEGM(xbaseline,Fbaseline,toPlot)
    if nargin < 3
       toPlot = 0;
    end
%     
% xbaseline = r{i}.x;
% Fbaseline= r{i}.y;
% toPlot = 1;



    xbaseline = xbaseline(find(Fbaseline>-50));
    Fbaseline = Fbaseline(find(Fbaseline>-50));


    y=Fbaseline;
    u = ones(size(y));
    for R2 = 5:5:15
        segm = segment([y,u],[0 1 1],R2);
        disp(R2)
        segmNotNan = segm(~isnan(segm));
        disp(var(segmNotNan(end-100:end)))
        if sum(isnan(segm)) < 500 % How much filtering to do?
            break
        end
    end
%     plot(segm)
    
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

    % Prune peaks that are not close to any point
    rmsLevel = 3*std(Fbaseline(end-100:end));
    [pks,locs] = findpeaks(newY,newX,'MinPeakProminence',5,'MinPeakDistance',8,'MinPeakHeight',12);
    newPks = [];
    newLocs = [];
    for i=1:length(pks)
%         disp(sprintf('(%2.0f,%2.0f)',locs(i),pks(i)))
        closestDistance = 1000;
        for j=1:length(xbaseline)
            theDist = pdist([[xbaseline(j),Fbaseline(j)];[locs(i),pks(i)]],'euclidean');
            if theDist < closestDistance
                closestDistance = theDist;
            end
        end
        if closestDistance < 50
            newPks = [newPks; pks(i)];
            newLocs = [newLocs; locs(i)];
        end
    end
    pks = newPks;
    locs = newLocs;
    
    % Prune peaks that are close to the start
    pks = pks(find(locs>5));
    locs = locs(find(locs>5));
    
    if toPlot > 0
        subplot(2,1,1)
        plot(xbaseline,Fbaseline,locs,pks,'or')
        subplot(2,1,2)
        plot(locs,pks,'or',newX,newY,'m-')
    end


end