clear all;
close all;
clear all;


function [matching] = getMatches(vec1,vec2)

end


allLs = []
vecnum = 1

fileList = getAllFiles('./');
for i=1:length(fileList)
    hasMat = findstr(char(fileList(i)),'s.mat');
    if length(hasMat)>0
       disp(sprintf('%s',char(fileList(i)))) 
       load(char(fileList(i)))
       rec = experiment.recording;
       for j=1:length(rec)
           if length(rec{j}.L)>1 && rec{j}.L(end)>130
%                 plot(rec{j}.x-rec{j}.L(end),rec{j}.y)
%                 hold on;
            t = rec{j}.L-mean(rec{j}.L(end-1:end))+160;
            vec{vecnum} = t;
            fvec{vecnum} = rec{j}.yPeaks;
            vecnum = vecnum + 1;
           end
       end
    end
end



distMatrix = zeros(length(vec),length(vec));
for ii=1:length(vec)
    for jj=1:length(vec)
        vec1m = flipud(sortrows([vec{ii} fvec{ii}]));
        vec2m = flipud(sortrows([vec{jj} fvec{jj}]));
        vec1 = sort(vec{ii},'descend');
        vec2 = sort(vec{jj},'descend');
        vec1 = reshape(vec1,1,length(vec1));
        vec2 = reshape(vec2,1,length(vec2));

        while length(vec1) < length(vec2)
            vec1 = [vec1 50];
            vec1m = [vec1m; 0 0];
        end
        while length(vec2) < length(vec1)
            vec2 = [vec2 50];
            vec2m = [vec2m; 0 0];
        end
        
        m = zeros(length(vec1),length(vec2));
        for n=1:length(vec1)
            for p=1:length(vec2)
                m(n,p) = abs(vec1(n)-vec2(p));
            end
        end

        totalError = 0;
        colsDone = zeros(length(vec1),1);
        rowsDone = zeros(length(vec2),1);

        while sum(colsDone)<(length(colsDone))*(length(colsDone)+1)/2
            mn = min(min(m));
            [rows,cols,vals] = find(m==mn);
            if length(rows)>1
                rows=rows(1);
                cols=cols(1);
            end
            if colsDone(cols)==0 && rowsDone(rows)==0
                colsDone(cols)=rows;
                rowsDone(rows)=cols;
                totalError = totalError + mn + 0*abs(vec1m(cols,2)-vec2m(rows,2));
            end
            m(rows,cols) = 500;
        end
        
        vec1mean = 0;
        vec2mean = 0;
        numVals = 0;
        for i=1:length(colsDone)
            if (rowsDone(i) > 0 && colsDone(i) > 0)
                vec1mean = vec1mean + vec1(rowsDone(i));
                vec2mean = vec2mean + vec2(colsDone(i));
                numVals = numVals + 1;
            end
        end
        vec1 = vec1 - vec1mean/numVals + 100;
        vec2 = vec2 - vec2mean/numVals + 100;
        
        
            
        m = zeros(length(vec1),length(vec2));
        for n=1:length(vec1)
            for p=1:length(vec2)
                m(n,p) = abs(vec1(n)-vec2(p));
            end
        end

        totalError = 0;
        colsDone = zeros(length(vec1),1);
        rowsDone = zeros(length(vec2),1);

        while sum(colsDone)<(length(colsDone))*(length(colsDone)+1)/2
            mn = min(min(m));
            [rows,cols,vals] = find(m==mn);
            if length(rows)>1
                rows=rows(1);
                cols=cols(1);
            end
            if colsDone(cols)==0 && rowsDone(rows)==0
                colsDone(cols)=rows;
                rowsDone(rows)=cols;
                totalError = totalError + mn +  0*abs(vec1m(cols,2)-vec2m(rows,2));
            end
            m(rows,cols) = 500;
        end    
        
        distMatrix(ii,jj) = totalError;
        distMatrxi(jj,ii) = totalError;
    end
        
end


Z = linkage(distMatrix,'ward','euclidean')


close all;
figure;
subplot(1,2,2)
dendrogram(Z,'Orientation','right')
cidx = cluster(Z,'MaxClust',9)




allBlocks = zeros(length(vec),250);
for i=1:max(cidx)
    vecmean{i} = 0
    hst{i} = []
end

for i=1:length(cidx)
%     vecmean{cidx(i)} = 100;
    vecmean{cidx(i)} = vecmean{cidx(i)} + mean(vec{i})/sum(cidx==cidx(i));
end

minVal = 0
for i=1:length(vec)
    pos{i} = vec{i}-vecmean{cidx(i)};
    for j=1:length(pos{i})
        foo = round(pos{i}(j))-1;
        if foo < minVal
            minVal = foo
        end
    end
end
for i=1:length(vec)
    pos{i} = vec{i}-vecmean{cidx(i)};
    for j=1:length(pos{i})
        hst{cidx(i)} = [hst{cidx(i)}; pos{i}(j)];
        pos{i}(j) = round(pos{i}(j))-minVal+2;
        allBlocks(i,pos{i}(j)) = 1;
        allBlocks(i,pos{i}(j)+1) = 1;
        allBlocks(i,pos{i}(j)-1) = 1;
    end
end



newBlock = zeros(size(allBlocks));
row = 1;
for i=max(cidx):-1:1
    for j=1:length(cidx)
        if cidx(j)==i
            newBlock(row,:) = allBlocks(j,:)*cidx(j);
            row = row + 1;
        end
    end
end
subplot(1,2,1)
imagesc(imrotate(newBlock,0))
xlabel('Contour length [nm]')
ylabel('Record num')


maxHst = 0;
maxLength = 0;
for i=1:length(hst)
    if length(hst{i}) > maxLength
        maxHst = i;
        maxLength = length(hst{i});
    end
end
figure;
hist(hst{maxHst},20)



vecnum = 1

fileList = getAllFiles('./');
for i=1:length(fileList)
    hasMat = findstr(char(fileList(i)),'s.mat');
    if length(hasMat)>0
       disp(sprintf('%s',char(fileList(i)))) 
       load(char(fileList(i)))
       rec = experiment.recording;
       for j=1:length(rec)
           if length(rec{j}.L)>1 && rec{j}.L(end)>130
               figure(cidx(vecnum)+10)
               pos2 = mean(rec{j}.L - vec{vecnum});
               pos2 = pos2 + mean(vec{vecnum}) - vecmean{cidx(vecnum)};
               x2 = rec{j}.x - pos2;
             plot(x2,rec{j}.y)
             hold on;
            vecnum = vecnum + 1;
           end
       end
    end
end
