close all;
clear all;

addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
addpath('../peaks/');
addpath('../benchmarks');


%% Use a simulated "real" data set
for i=1:20
    % Full Luciferase with random i27
    numI27 = 4%randi(9,1,1)-1;
    params.Force = sort([60 60 60 200*ones(1,numI27)]);
    params.ForceSD = 7.5;
    params.Lc = sort([30 80 130 190:28:190+28*(numI27-1)]);
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
    
    [pks,locs]=getPeaksLC(r{i}.x,r{i}.y,0);
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
for i=21:40
    % Full Luciferase split with random i27
    numI27 =4%randi(9,1,1)-1;
    params.Force = sort([60 60 60 60 200*ones(1,numI27) ]);
    params.ForceSD = 7.5;
    params.Lc = sort([30 80 100 130 190:28:190+28*(numI27-1) ]);
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
% for i=16:20
%     % Full Luciferase with random i27
%     numI27 =randi(9,1,1)-1;
%     params.Force = [20 30 40 200*ones(1,numI27)];
%     params.ForceSD = 7.5;
%     params.Lc = sort([30 80 130 190:28:190+28*(numI27-1)]);
%     params.LcSD = 1;
%     params.Persistence = 0.4;
%     params.rmsNoise = 10;
%     params.nonSpecificForce = 50;
%     params.nonSpecificForceSD = 10;
%     params.nonSpecificForces = 1;
%     params.curveLength = 450;
%     params.diffusion = 0.0006;
%     [data, Foriginal,rlocs,rpeaks] = generateSpecial(params);
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
originalR = r;


r=originalR;
[r distMatrix] = calculatingMatrix2(r);


figure(10)
subplot(2,1,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment22(r,ones(size(r)),1,10,ones(size(r)));
ylabel('Record Number')

subplot(2,1,2)
num = 0
Ldist = [];
for i=1:length(r)
    for j=1:length(r{i}.L)
        Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR;
        Ldist = [Ldist; Ldistcalc];
    end
end
[f,xi]=ksdensity(Ldist,0:1:max(Ldist),'bandwidth',4);
f = f/max(f)
num = num + 1;
plot(xi,f)
axis([0 max(Ldist) 0 max(f)*1.2])
xlabel('Contour length')

for i=1:length(r)
    r{i}.consistency = zeros(size(r{i}.L));
    for j=1:length(r{i}.L)
        Ldistcalc = (r{i}.L(j) + r{i}.Ladj)-minR;
        tmp = abs(xi-Ldistcalc);
        [idx idx] = min(tmp); %index of closest value
        closest = xi(idx); %closest value
        r{i}.consistency(j) = f(closest);
    end
end    


[r distMatrix] = calculatingMatrix3(r);
[r distMatrix] = calculatingMatrix(r);
[r distMatrix] = calculatingMatrix2(r);
[r distMatrix] = calculatingMatrix2Dlc(r);
figure(11)
numClusters = 3;
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
group = cidx;
subplot(1,2,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment22(r,cidx,perm,10,group);
disp(numClusters)
newCriterion = log(mean(meanDifferenceSDmean)) + log(numClusters);
ylabel('Record Number')
xlabel('Contour length')
