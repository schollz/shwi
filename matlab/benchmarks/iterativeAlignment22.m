function [r,newMinR, meanDifferenceSDmean] = iterativeAlignment22(r,cidx,perm,iterations,group)
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
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (30--10)*rand(1,1);
%     cidx(i) = 1;
% end
% for i=21:40
%     r{i}.L = [ 30 + lcSD.*randn(1,1) 70 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (30--10)*rand(1,1);
%     cidx(i) = 1;
% end
% % for i=41:60
% %     r{i}.L = [ 30 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)] + -10 + (30--10)*rand(1,1);
% %     cidx(i) = 1;
% % end
% perm = 1:max(cidx);
% iterations = 1;
% group = cidx;


    % Initialize matching
    tic
    textprogressbar('matching peaks: ')
    num = 0
    for i=1:length(r)
        for j=1:length(r)
            num = num + 1;
            textprogressbar(num/(length(r)^2)*100)
            if i~=j
                r{i}.bestMatch{j} = getBestMatch(r{i}.L,r{j}.L);
            end
        end
    end
    textprogressbar('done.')

    % Initialize Ladj for each
    for i=1:length(r)
        r{i}.Ladj = 0;
    end

    % Iterative alignment
    drawPieces(r)
    pause(0.1)
%    textprogressbar('iterative alignment: ');
    isDone = zeros(length(cidx),1); 
    for ci=1:max(cidx)
        lastMovementMean = 1000;
        lastMovement = [0 0 0];
        movementType=ones(1000,1)*1;
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
                        [shiftAmount] = getAlignmentDifference2(r{j}.L-r{j}.Ladj,r{i}.L-r{i}.Ladj,r{j}.bestMatch{i});
                        newAdj = [newAdj; shiftAmount];  
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
            disp(sprintf('%2.3f %2.3f\n',mean(lastMovement),abs(mean(lastMovement) - lastMovementMean)))
            if abs(mean(lastMovement) - lastMovementMean) < 1
                break
            end
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





