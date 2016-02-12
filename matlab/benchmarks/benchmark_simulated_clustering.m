addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
close all;
clear all;
warning('off','all')
warning
params.Persistence = 0.4;
% 
% %% Use a simulated data set
% for i=1:50
%     params.Force = 100;
%     params.ForceSD = 1;
%     params.Lc = 28;
%     params.LcSD = 3;
%     params.numPeaks = randi([7 12],1,1);
%     params.rmsNoise = 10;
%     params.nonSpecificForce = 50;
%     params.nonSpecificForceSD = 10;
%     params.nonSpecificForces = 0;
%     params.curveLength = 450;
%     params.diffusion = 0.0003;
%     [data, Foriginal,rlocs,rpeaks] = generateCurve(params);
%     
%     [origin] = getOrigin(data(:,1),data(:,2),0);
%     r{i}.x = data(:,1)-origin(:,1);
%     r{i}.y = data(:,2)-origin(:,2);
%     clear data
%     r{i}.file = sprintf('%d',i);
%     r{i}.name = sprintf('simulated');
%     r{i}.realL = rlocs;
%     r{i}.realF = rpeaks;
%     
%     [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,0);
%     r{i}.xPeaks = locs;
%     r{i}.F = pks;
%     r{i}.L=[];
%     peaks = [];
%     for k=1:length(r{i}.xPeaks)
%         x=r{i}.xPeaks(k);
%         F=r{i}.F(k);
%         L0 = getLc(params.Persistence,x,F);
%         r{i}.L = [r{i}.L; L0];
%     end
%     pause(0.1)
%     plot(r{i}.x,r{i}.y,rlocs,rpeaks,'x',locs,pks,'o')
% end
% originalR = r;


%% Use a real data
filelist1 = getAllFiles('D:\Marszalek Lab\Force Curves\AllYeastPGK\YeastPGK');
filelist1 = getAllFiles('../testdata/');
i = 0;
clear filelist
clear r
for num=1:length(filelist1)
    if ~isempty(findstr(filelist1{num},'.afm'))==1
        i = i + 1;
        try
            data = importdata(char(filelist1{num}));
            [origin] = getOrigin(data(:,1),data(:,2),0);
            r{i}.x = data(:,1)-origin(:,1);
            r{i}.y = data(:,2)-origin(:,2);
            clear data
            r{i}.file = sprintf('%d',i);
            r{i}.name = char(filelist1{num});

            [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,0);
            r{i}.xPeaks = locs;
            r{i}.F = [];
            r{i}.L=[];
            peaks = [];
            for k=1:length(r{i}.xPeaks)
                x=r{i}.xPeaks(k);
                F=pks(k);
                try
                    L0 = getLc(params.Persistence,x,F);
                    if ~isnan(L0)
                        r{i}.L = [r{i}.L; L0];
                        r{i}.F = [r{i}.F; pks(k)];
                    end
                catch
                end
            end
            pause(0.1)
            plot(r{i}.x,r{i}.y,locs,pks,'o')
            if length(r{i}.L) < 2
                i = i - 1;
            end
        catch
            i = i - 1;
        end
    end
end
originalR = r;

r = originalR;

%% Iterative pruning
% Calculate distance matrix
[r distMatrix] = calculatingMatrix(r);

% Align 
subplot(2,1,1)
numClusters = 2
Z = linkage(distMatrix,'ward','euclidean');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
group = cidx;
subplot(2,1,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,200,group);

% Calculate Kernel density
for k=1:max(cidx)
    Ldist = []
    for i=1:length(r)
        for j=1:length(r{i}.L)
            if cidx(i)==k
                Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR+10;
                Ldist = [Ldist; Ldistcalc];
            end
        end
    end
    [f,xi]=ksdensity(Ldist,[sort(unique(Ldist))]);
    thresholdDensity = mean(f(find(xi>mean(Ldist)+2*std(Ldist))));
    subplot(2,1,2)
    plot(xi,f)
    axis([0 470 0 6e-3])
    for i=1:length(r)
        if cidx(i)==k
            newRL = [];
            newXpeaks = [];
            newF = [];
            for j=1:length(r{i}.L)
                Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR+10;
                if f(xi==Ldistcalc) > thresholdDensity
                    newRL =[newRL; r{i}.L(j)];
                    newXpeaks = [newXpeaks; r{i}.xPeaks(j)];
                    newF = [newF; r{i}.F(j)];
                else
                    disp(r{i}.L(j)-r{i}.Ladj)
                end
            end
            r{i}.F = newF;
            r{i}.xPeaks = newXpeaks;
            r{i}.L = newRL;
        end
    end
end

for i=1:length(r)
    if length(r{i}.L)==0
        disp(i)
    end
end




%% AIC to find number of clusters
[r distMatrix] = calculatingMatrix(r);
    
% AIC criterion to find number of clusters
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
lastCriterion = 10000;
for numClusters=2:1:15
    [H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
    ylabel('Group Number')
    cidx = cluster(Z,'MaxClust',numClusters);
    cidx = T;
    group = cidx;
    subplot(1,2,1)
    [r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,30,group);
    disp(numClusters)
    newCriterion = -2*log(1/mean(meanDifferenceSDmean)) + 2*numClusters;  %AIC criterion
    disp(newCriterion)
    if newCriterion > lastCriterion
%         break
    end
    lastCriterion = newCriterion;
    ylabel('Record Number')
    xlabel('Contour length')
end
numClusters = numClusters -1


%% Actually plot clusters
numClusters = 7
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
group = cidx;
subplot(1,2,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,200,group);
disp(numClusters)
newCriterion = log(mean(meanDifferenceSDmean)) + log(numClusters);
ylabel('Record Number')
xlabel('Contour length')





% Plot one
ci=6;
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
