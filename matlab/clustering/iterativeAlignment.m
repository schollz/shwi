function [r,newMinR, meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,iterations,group)
    %% Perform iterative alignment
clear r
clear cidx 
clear perm
clear group
trueCidx = [];
lcSD=1;

% for i=1:20
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (80--10)*rand(1,1);
%     cidx(i) = 1;
% end
% for i=20:50
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1) 100.*randn(1,1)] + -10 + (80--10)*rand(1,1);
%     cidx(i) = 1;
% end
% for i=50:60
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (80--10)*rand(1,1);
%     cidx(i) = 1;
% end

for i=1:25
    r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (80--10)*rand(1,1);
    cidx(i) = 1;
end
for i=26:50
    r{i}.L = [ 12 + lcSD.*randn(1,1) 24 + lcSD.*randn(1,1) 36 + lcSD.*randn(1,1) 100.*randn(1,1)] + -10 + (80--10)*rand(1,1);
    cidx(i) = 1;
end

perm = 1:max(cidx);
iterations = 1;
group = cidx;
perm = 1:max(cidx);
iterations = 30;
group = cidx;
[r distMatrix] = calculatingMatrix(r);

    % Initialize Ladj for each
    for i=1:length(r)
        r{i}.Ladj = 0;
    end
    drawPiecesGetIteration1(r)
    pause(0.1)
    % Iterative alignment
    tic
    textprogressbar('iterative alignment: ');
    numIterations = iterations;
    for numRuns = 1:numIterations
        textprogressbar((numRuns/numIterations*100))
        for ci=1:max(cidx)
            cind = find(cidx==ci);
            meanDifferenceSD{ci} = [];
            for jj=1:length(cind)
                i = cind(jj);
                meanDifferences = zeros(length(cind),1);
                for jj2=1:length(cind)
                    if (jj==jj2)

                    else
                        i2 = cind(jj2);
                        if sum(r{i}.L) == sum(r{i2}.L)
                            meanDifferences(jj2) = 0;
                        else
                            meanDifferences(jj2)=-1*getDifference(r{i}.L+r{i}.Ladj,r{i2}.L+r{i2}.Ladj,r{i}.matches{i2});
                        end
                    end
                end
                meanDifference = mean(meanDifferences(meanDifferences>0));
                meanDifferenceSD{ci} = [meanDifferenceSD{ci}; std(meanDifferences)];
                if isnan(meanDifference)==0
                    r{i}.Ladj = r{i}.Ladj + meanDifference;
                end

            end
        end
            drawPiecesGetIteration1(r)
            pause(0.05)
    end
    textprogressbar('done');
    toc
    
    meanDifferenceSDmean = zeros(max(cidx),1);
    for i=1:length(meanDifferenceSD)
        if sum(meanDifferenceSD{i})==0 || length(meanDifferenceSD{i})<0.1*length(cidx)
            meanDifferenceSDmean(i)=0;
        else
            meanDifferenceSDmean(i) = mean(meanDifferenceSD{i});
        end
    end
    for i=1:length(meanDifferenceSDmean)
        if meanDifferenceSDmean(i)==0
            meanDifferenceSDmean(i) = max(meanDifferenceSDmean);
        end
    end

    % Align to last peak for each
    meanAll = zeros(max(cidx),1);
    for ci=1:max(cidx)
        cind = find(cidx==ci);
        for jj=1:length(cind)
            i = cind(jj);
            meanAll(ci) = meanAll(ci) + (r{i}.L(end)+r{i}.Ladj(end))/length(cind);      
        end
    end

    minR = 1000;
    maxR = -1000;
    for ci=1:max(cidx)
        cind = find(cidx==ci);
        for jj=1:length(cind)
            i = cind(jj);
            r{i}.Ladj = r{i}.Ladj - meanAll(ci);
            if max((r{i}.L + r{i}.Ladj)) > maxR
                maxR = max((r{i}.L + r{i}.Ladj));
            end
            if min((r{i}.L + r{i}.Ladj)) < minR
                minR = min((r{i}.L + r{i}.Ladj));
            end      
        end
        
    end



    allBlocks = zeros(length(r)+length(cind)+1,round(maxR)-round(minR)+10);
    allBlocks = zeros(length(r),300);
    row = 1;
    for ci2=1:max(cidx)
        ci = perm(ci2);
        cind = find(cidx==ci);
        for jj=1:length(cind)
            i = cind(jj);
            for j=1:length(r{i}.L)
                val = round((r{i}.L(j) + r{i}.Ladj))-round(minR)+10;
%                 allBlocks(row,val) = ci;
%                 allBlocks(row,val+1) = ci;
%                 allBlocks(row,val-1) = ci;  
                allBlocks(row,val) = group(i);
                allBlocks(row,val+1) = group(i);
                allBlocks(row,val-1) = group(i);  
            end
            row = row + 1;
        end
        allBlocks(row,:) = -99;
        row = row + 1;
    end
    colormap jet
    allBlocks(allBlocks==0) = 1+max(max(allBlocks));
    allBlocks(allBlocks==-99) = 1+max(max(allBlocks));
    imagesc(flipud(allBlocks))
    newMinR = -1 * (-round(minR)+10);
end