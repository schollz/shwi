function [r,newMinR, meanDifferenceSDmean] = iterativeAlignment2(r,cidx,perm,iterations,group)
%% Perform iterative alignment
% 
% 
% clear r
% clear cidx 
% clear perm
% clear group
% trueCidx = [];
% lcSD=1;
% 
% for i=1:20
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1) 5*rand(1,1)] + -10 + (30--10)*rand(1,1);
%     cidx(i) = 1;
% end
% for i=1:20
%     r{i}.L = [ 12 + lcSD.*randn(1,1) 24 + lcSD.*randn(1,1) 36 + lcSD.*randn(1,1) 40.*randn(1,1)] + -10 + (80--10)*rand(1,1);
%     cidx(i) = 1;
% end
% perm = 1:max(cidx);
% iterations = 1;
% group = cidx;



    % Initialize Ladj for each
    for i=1:length(r)
        r{i}.Ladj = 0;
    end
tic
    % Iterative alignment
    drawPieces(r)
    pause(0.1)
%    textprogressbar('iterative alignment: ');
    isDone = zeros(length(cidx),1); 
    for ci=1:max(cidx)
        lastMovementMean = 1000;
        lastMovement = [0 0 0];
        movementType=[3 3 3 3 3 3 2 2 2 2 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
        for kkk=1:length(movementType)
%             textprogressbar((kkk+20*(ci-1))/(20*max(cidx))*100)
            movement = 0;
            toProcess =find(cidx==ci);
            iArray = toProcess(randperm(length(toProcess)));
            jArray = toProcess(randperm(length(toProcess)));
            for ii=1:length(iArray)-1
                i=iArray(ii);
                newAdj = [];
                for jj=ii:length(jArray)
%                     disp(sprintf('%d %d %d %d',ci,kkk,ii,jj))
                    j=jArray(jj);
                    if j~=i
                        try
                        [shiftAmount] = getAlignmentDifference(r{j}.L-r{j}.Ladj,r{i}.L-r{i}.Ladj);
                        newAdj = [newAdj; shiftAmount];  
                        catch
                        end
                    end
                end
                if mean(newAdj) > 0
                    r{i}.Ladj = r{i}.Ladj + movementType(kkk);
                    movement = movement + mean(newAdj);
                elseif mean(newAdj) < 0
                    r{i}.Ladj = r{i}.Ladj - movementType(kkk);
                    movement = movement - mean(newAdj);
                end
            end
            drawPieces(r)
            pause(0.1)
            lastMovement = [movement lastMovement];
            lastMovement = lastMovement(1:3);
            if abs(mean(lastMovement) - lastMovementMean) < 1
                break
            end
            disp(sprintf('%2.0f %2.0f\n',mean(lastMovement),abs(mean(lastMovement) - lastMovementMean)))
            lastMovementMean = mean(lastMovement);
        end
    end
%     textprogressbar('done');
    toc
    for i=1:length(r)
        r{i}.Ladj = -1*r{i}.Ladj;
    end
    
    
    minR = 1000;
    maxR = -1000;
    for ci=1:max(cidx)
        cind = find(cidx==ci);
        for jj=1:length(cind)
            i = cind(jj);
            r{i}.Ladj = r{i}.Ladj;
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
    meanDifferenceSDmean = 0;
end





