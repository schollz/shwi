addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
close all;
clear all;
warning('off','all')
warning


folder = '../testdata/8i27';
r = collectMats(folder,0);


% analyzeInitialLengths(r)

% Simulate contour length increments
% clear r
trueCidx = []
for i=1:20
    r{i}.L = [ 10 + 4.*randn(1,1) 40 + 4.*randn(1,1) 100 + 4.*randn(1,1)]
    trueCidx(i) = 1;
end
for i=20:40
    r{i}.L = [ 10 + 4.*randn(1,1) 40 + 4.*randn(1,1) 100 + 4.*randn(1,1)]
    trueCidx(i) = 1;
end
for i=40:60
    r{i}.L = [ 10 + 4.*randn(1,1) 40 + 4.*randn(1,1) 60 + 4.*randn(1,1) 100 + 4.*randn(1,1)]
    trueCidx(i) = 3;
end
for i=60:70
    r{i}.L = [ 10 + 7.*randn(1,1)  30 + 7.*randn(1,1) 170 + 7.*randn(1,1)]
end

%% Calculate distance matrix
[r distMatrix] = calculatingMatrix2(r);

%% Plot mdscaled curves
% Y = mdscale(distMatrix,2)
% figure(123)
% for i=1:size(Y,1)
%     x = r{i}.x;
%     y = r{i}.y;
%     ind = find(y>-10);
%     x = x(ind);
%     y = y(ind);
%     x = x - mean(x);
%     y = y-min(y);
%     x = x/20;
%     y = y/20;
%     hold on;
%     plot(x+Y(i,1),y+Y(i,2))
% end

%% Grouping based on PGK
% group = zeros(length(r),1);
% for i=1:length(r)
%     if ~isempty(findstr(char(r{i}.name),'YeastPGK'))
%         group(i) = 1;
%     elseif  ~isempty(findstr(char(r{i}.name),'yPGKccCT'))
%         group(i) = 3;
%     elseif  ~isempty(findstr(char(r{i}.name),'yPGKCCdumped'))
%         group(i)=2;
%     elseif  ~isempty(findstr(char(r{i}.name),'yPGK68cc'))
%         group(i)=4;
%     elseif  ~isempty(findstr(char(r{i}.name),'PGK88'))
%         group(i)=5;
%     elseif  ~isempty(findstr(char(r{i}.name),'CysYeastPGK'))
%         group(i)=1;
%     end
% end


%% Determine number of clusters
% AIC criterion to find number of clusters
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
lastCriterion = 10000;
for numClusters=11:1:15
    [H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
    ylabel('Group Number')
    cidx = cluster(Z,'MaxClust',numClusters);
    cidx = T;
    group = cidx;
    subplot(1,2,1)
    [r,minR,meanDifferenceSDmean] = iterativeAlignment22(r,cidx,perm,30,group);
    disp(numClusters)
    newCriterion = -2*log(1/mean(meanDifferenceSDmean)) + 2*numClusters;  %AIC criterion
    disp(newCriterion)
    if newCriterion > lastCriterion
        break
    end
    lastCriterion = newCriterion;
    ylabel('Record Number')
    xlabel('Contour length')
end
numClusters = numClusters -1


Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
group = cidx;
subplot(1,2,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,100,group);
disp(numClusters)
newCriterion = log(mean(meanDifferenceSDmean)) + log(numClusters);
ylabel('Record Number')
xlabel('Contour length')
    

%% Analyze efficiency of assignment (for in 
for i=1:max(trueCidx)
    subCidx = cidx(find(trueCidx==i));
    percentWrong = length(find(subCidx ~= mode(subCidx)))/length(subCidx) 
end
    




% Plot one
ci=1;
cind = find(cidx==ci);
cm=colormap(jet(length(cind)+1));
num = 0
for jj=1:length(cind)
    num = num + 1;
    i = cind(jj);
    figure(10)
    plot(r{i}.x+r{i}.Ladj,r{i}.y,'color',cm(num,:))
    axis([-350 250 -200 400])
    hold on;
end






% D:\Marszalek Lab\Force Curves\yPGKCCdumped
forces = [];
L = 32+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{1} = forces;
forces = [];
L = 104+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{2} = forces;
forces = [];
L = 155+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 15 && Ls(j) < L + 15)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{3} = forces;
forces = [];
L = 65+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 9 && Ls(j) < L + 9)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{4} = forces;
save('../PGK-analysis/yPGKCCdumped.mat','yPGKCCdumped')


% D:\Marszalek Lab\Force Curves\PGK88
forces = [];
lengths=[];
L = 30+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
           lengths = [lengths; Ls(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
forces = [];
L = 72+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 15 && Ls(j) < L + 15)
           forces = [forces; r{i}.F(j)];
           lengths = [lengths; Ls(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{1} = forces;
forces = [];
L = 114+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{2} = forces;
forces = [];
L = 160+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKCCdumped{3} = forces;
save('../PGK-analysis/PGK88.mat','yPGKCCdumped')




% D:\Marszalek Lab\Force Curves\yPGKccCT
forces = [];
L = 45+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKccCT{1} = forces;
forces = [];
L = 122+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKccCT{2} = forces;
forces = [];
L = 174+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKccCT{3} = forces;
forces = [];
L = 88+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGKccCT{4} = forces;
save('../PGK-analysis/yPGKccCT.mat','yPGKccCT')



% D:\Marszalek Lab\Force Curves\YeastPGK
forces = [];
L = 22+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 15 && Ls(j) < L + 15)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
YeastPGK{1} = forces;
forces = [];
L = 74+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
YeastPGK{2} = forces;
forces = [];
L = 125+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
YeastPGK{3} = forces;
save('../PGK-analysis/YeastPGK.mat','YeastPGK')



% D:\Marszalek Lab\Force Curves\yPGK68cc
forces = [];
L = 55+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 15 && Ls(j) < L + 15)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGK68cc{1} = forces;
forces = [];
L = 106+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGK68cc{2} = forces;
forces = [];
L = 150+minR;
for i=1:length(r)
    Ls = r{i}.L + r{i}.Ladj;
    for j=1:length(Ls)
       if (Ls(j) > L - 10 && Ls(j) < L + 10)
           forces = [forces; r{i}.F(j)];
       end
    end
end
disp(sprintf('Peak %d: %2.1f +/- %2.1f',L,mean(forces),std(forces)))
yPGK68cc{3} = forces;
save('../PGK-analysis/yPGK68cc.mat','yPGK68cc')










clear peaks;
groups = [6 13];
numPeaks = 4;
for j=1:numPeaks
    peaks{j} = [];
end
for i=1:length(r)
    if (any(abs(group(i)-groups)<1e-10) && length(r{i}.F)>=numPeaks+1)
        for j=1:numPeaks
            peaks{j} = [peaks{j}; r{i}.F(j)];
        end
    end
end
for j=1:numPeaks
    disp(sprintf('Peak %d: %2.1f +/- %2.1f',j,mean(peaks{j}),std(peaks{j})))
end
figure(21)
nhist(peaks)




%% Plot all
for ci=1:max(cidx)
    cind = find(cidx==ci);
    num = 0
    cm=colormap(jet(1+length(cind)));
    for jj=1:length(cind)
        num = num + 1;
        i = cind(jj);
        figure(10+1)
        plot(r{i}.x+r{i}.Ladj,r{i}.y,'color',cm(num,:))
        hold on;
    end
end


% Plot one
ci=1;
cind = find(cidx==ci);
cm=colormap(jet(length(cind)+1));
num = 0
for jj=1:length(cind)
    num = num + 1;
    i = cind(jj);
    figure(10)
    plot(r{i}.x+r{i}.Ladj,r{i}.y,'color',cm(num,:))
    axis([-350 250 -200 400])
    hold on;
end

% Plot two
cind = find(cidx==1 | cidx==3);
cm=colormap(jet(length(cind)+1));
num = 0
for jj=1:length(cind)
    num = num + 1;
    i = cind(jj);
    figure(10)
    plot(r{i}.x+r{i}.Ladj,r{i}.y,'color',cm(num,:))
    disp([r{i}.name,r{i}.file])
    axis([-350 250 -200 400])
    hold on;
end

% 
% % %% Plot one Lcs/F
% num = 0;
% cm=colormap(jet(numClusters+1));
% for ci=1:numClusters
% cind = find(cidx==ci);
%     num = num + 1;
%     figure(400);
% for jj=1:length(cind)
%     i = cind(jj);
%     plot((r{i}.L+r{i}.Ladj),r{i}.F(1:end),'o-','color',cm(num,:))
%     hold on;
%     axis([-200 20 0 100])
% end
% end
% 
% % Hist Some
% allLs = [];
% goodGroups = [1:numClusters]
% for ci=1:length(goodGroups)
%     cind = find(cidx==goodGroups(ci))
%     for jj=1:length(cind)
%         i = cind(jj);
%         for kk=length(r{i}.L)-1
%             allLs = [allLs; r{i}.L(kk)+r{i}.Ladj];
%         end
%     end
% end
% figure(571)
% hist(allLs,10)
% load ../a.mat
% a{1}=allLs;
% save('../a.mat','a')
% 
% %% Hist all
% allAllLs = []
% for ci=1:max(cidx)
%     cind = find(cidx==ci);
%     allLs = [];
%     numPeaks = [];
%     for jj=1:length(cind)
%         i = cind(jj);
%         for kk=1:length(r{i}.L)
%             allLs = [allLs; r{i}.L(kk)+r{i}.Ladj];
%             allAllLs = [allAllLs;  r{i}.L(kk)+r{i}.Ladj];
%         end
%         numPeaks = [numPeaks; length(r{i}.L)];
%     end
%     infrData{ci}.Ls = allLs;
%     infrData{ci}.peaks = mean(numPeaks);    
% end
% for i=1:length(infrData)
%     subplot(length(infrData),1,i)
%     hist(infrData{i}.Ls,50)
%     [mu,sig,p]=gaussian_mixture_model(infrData{i}.Ls',round(infrData{i}.peaks),.0001)
%     diff(sort(mu))
% end
% figure;
% hist(Ls{1})
% 
% 
% 
% r = iterativeAlignment(r,1:length(cidx),1:length(cidx),1);
% ylabel('Record Number')
% xlabel('Contour length')
% 
