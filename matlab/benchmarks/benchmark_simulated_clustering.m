addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
addpath('../peaks/');
close all;
clear all;
warning('off','all')
warning
params.Persistence = 0.4;

%% Use a simulated data set
for i=1:10
    params.Force = 100;
    params.ForceSD = 1;
    params.Lc = 28;
    params.LcSD = 3;
    params.numPeaks = randi([7 12],1,1);
    params.rmsNoise = 5;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 0;
    params.curveLength = 450;
    params.diffusion = 0.0003;
    [data, Foriginal,rlocs,rpeaks] = generateCurve(params);
    
    [origin] = getOrigin(data(:,1),data(:,2),0);
    r{i}.x = data(:,1)-origin(:,1);
    r{i}.y = data(:,2)-origin(:,2);
    clear data
    r{i}.file = sprintf('%d',i);
    r{i}.name = sprintf('simulated');
    r{i}.realL = rlocs;
    r{i}.realF = rpeaks;
    
    [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,1);
    pause(0.1);
    r{i}.xPeaks = locs;
    r{i}.F = pks;
    r{i}.L=[];
    peaks = [];
    for k=1:length(r{i}.xPeaks)
        x=r{i}.xPeaks(k);
        F=r{i}.F(k);
        L0 = getLc(params.Persistence,x,F);
        r{i}.L = [r{i}.L; L0];
    end
end
originalR = r;


%% Use a simulated "real" data set
for i=1:5
    % Full Luciferase with Nonspecific
    params.Force = [20 30 40 200 200 200 200 200];
    params.ForceSD = 7.5;
    params.Lc = [30 80 130 190 220 250 280 310];
    params.LcSD = 1;
    params.Persistence = 0.4;
    params.rmsNoise = 10;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 1;
    params.curveLength = 450;
    params.diffusion = 0.0006;
    [data, Foriginal,rlocs,rpeaks] = generateSpecial(params);
    
    [origin] = getOrigin(data(:,1),data(:,2),0);
    r{i}.x = data(:,1)-origin(:,1);
    r{i}.y = data(:,2)-origin(:,2);
    clear data
    r{i}.file = sprintf('%d',i);
    r{i}.name = sprintf('simulated');
    r{i}.realL = rlocs;
    r{i}.realF = rpeaks;
    
    [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,0);
    r{i}.xPeaks = locs;
    r{i}.F = pks;
    r{i}.L=[];
    peaks = [];
    for k=1:length(r{i}.xPeaks)
        x=r{i}.xPeaks(k);
        F=r{i}.F(k);
        L0 = getLc(params.Persistence,x,F);
        r{i}.L = [r{i}.L; L0];
    end
    pause(0.1)
    plot(r{i}.x,r{i}.y,rlocs,rpeaks,'x',locs,pks,'o')
end
for i=6:20
    % Full Luciferase without Nonspecific
    params.Force = [20 30 40 200 200 200 200 200];
    params.ForceSD = 7.5;
    params.Lc = [30 80 130 190 220 250 280 310];
    params.LcSD = 1;
    params.Persistence = 0.4;
    params.rmsNoise = 10;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 0;
    params.curveLength = 450;
    params.diffusion = 0.0006;
    [data, Foriginal,rlocs,rpeaks] = generateSpecial(params);
    
    [origin] = getOrigin(data(:,1),data(:,2),0);
    r{i}.x = data(:,1)-origin(:,1);
    r{i}.y = data(:,2)-origin(:,2);
    clear data
    r{i}.file = sprintf('%d',i);
    r{i}.name = sprintf('simulated');
    r{i}.realL = rlocs;
    r{i}.realF = rpeaks;
    
    [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,0);
    r{i}.xPeaks = locs;
    r{i}.F = pks;
    r{i}.L=[];
    peaks = [];
    for k=1:length(r{i}.xPeaks)
        x=r{i}.xPeaks(k);
        F=r{i}.F(k);
        L0 = getLc(params.Persistence,x,F);
        r{i}.L = [r{i}.L; L0];
    end
    pause(0.1)
    plot(r{i}.x,r{i}.y,rlocs,rpeaks,'x',locs,pks,'o')
end
for i=21:30
    % Only last two peaks of Lucferase, without Nonspecific
    params.Force = [30 40 200 200 200 200 200];
    params.ForceSD = 5;
    params.Lc = [80 130 190 220 250 280 310];
    params.LcSD = 1;
    params.Persistence = 0.4;
    params.rmsNoise = 10;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 0;
    params.curveLength = 450;
    params.diffusion = 0.0006;
    [data, Foriginal,rlocs,rpeaks] = generateSpecial(params);
    
    [origin] = getOrigin(data(:,1),data(:,2),0);
    r{i}.x = data(:,1)-origin(:,1);
    r{i}.y = data(:,2)-origin(:,2);
    clear data
    r{i}.file = sprintf('%d',i);
    r{i}.name = sprintf('simulated');
    r{i}.realL = rlocs;
    r{i}.realF = rpeaks;
    
    [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,0);
    r{i}.xPeaks = locs;
    r{i}.F = pks;
    r{i}.L=[];
    peaks = [];
    for k=1:length(r{i}.xPeaks)
        x=r{i}.xPeaks(k);
        F=r{i}.F(k);
        L0 = getLc(params.Persistence,x,F);
        r{i}.L = [r{i}.L; L0];
    end
    pause(0.1)
    plot(r{i}.x,r{i}.y,rlocs,rpeaks,'x',locs,pks,'o')
end
originalR = r;



%% Use a real data
filelist1 = getAllFiles('D:\Marszalek Lab\Force Curves\3I27_Luciferase_4I27')
filelist1 = getAllFiles('D:\Marszalek Lab\Force Curves\AllYeastPGK\yPGKCCdumped');
filelist1 = getAllFiles('../testdata/');
filelist1 = getAllFiles('D:\Marszalek Lab\Force Curves\AllYeastPGK\YeastPGK');
i = 0;
clear filelist
clear r
for num=1:length(filelist1)
    if ~isempty(findstr(filelist1{num},'.afm'))==1
        i = i + 1;
        try
            disp(num)
            data = importdata(char(filelist1{num}));
            [origin] = getOrigin(data(:,1),data(:,2),1);
            pause(0.05)
            r{i}.x = data(:,1)-origin(:,1);
            r{i}.y = data(:,2)-origin(:,2);
            clear data
            r{i}.file = sprintf('%d',i);
            r{i}.name = char(filelist1{num});

            [pks,locs]=getPeaksSEGM(r{i}.x,r{i}.y,1,1);
            pause(0.6)
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
            if length(r{i}.L) < 2
                i = i - 1;
                disp(sprintf('skipping %d',i))  
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
[r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,1:max(cidx),50,group);
% [r,minR,meanDifferenceSDmean] = iterativeAlignment(r,ones(size(r)),ones(size(r)),200,ones(size(r)));

% Calculate Kernel density
for k=1:max(cidx)
    Ldist = [];
    for i=1:length(r)
        for j=1:length(r{i}.L)
            if cidx(i)==k
                Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR+10;
                Ldist = [Ldist; Ldistcalc];
            end
        end
    end
    [f,xi]=ksdensity(Ldist,[sort(unique(Ldist))],'bandwidth',5);
    thresholdDensity = mean(f)-std(f);
    if sum(cidx==k) < 5
        thresholdDensity = -1;
    end
    subplot(2,1,2)
    plot(xi,f)
    axis([0 700 0 20e-3])
    hist(Ldist,0:5:500)
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

% Remove any that were fully pruned
rOld = r;
clear r
new = 0;
for i=1:length(rOld)
    if length(rOld{i}.L)~=0 
        new = new + 1;
        r{new} = rOld{i};
    else
        disp(sprintf('removed %d',i))
    end
end





%% AIC to find number of clusters
[r distMatrix] = calculatingMatrix(r);
    
% AIC criterion to find number of clusters
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
lastCriterion = 10000;
for numClusters=2:1:6
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
figure(10)
numClusters = 4
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
group = cidx;
subplot(1,2,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment2(r,cidx,perm,100,group);
disp(numClusters)
newCriterion = log(mean(meanDifferenceSDmean)) + log(numClusters);
ylabel('Record Number')
xlabel('Contour length')

% Show densities
figure(12)
maxNum = 0
for k=1:max(cidx)
    if sum(cidx==k) > 5
        maxNum = maxNum + 1
    end
end
num = 0
for k=1:max(cidx)
    Ldist = [];
    for i=1:length(r)
        for j=1:length(r{i}.L)
            if cidx(i)==k
                Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR+10;
                Ldist = [Ldist; Ldistcalc];
            end
        end
    end
    [f,xi]=ksdensity(Ldist,0:1:max(Ldist),'bandwidth',4);
    if sum(cidx==k) > 3
        num = num + 1;
        subplot(maxNum,1,num)
        plot(xi,f)
        axis([0 max(Ldist) 0 max(f)*1.2])
    end
end




% Plot one
ci=3;
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
