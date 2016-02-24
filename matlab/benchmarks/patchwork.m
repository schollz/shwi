addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
addpath('../peaks/');
addpath('../benchmarks');

% params.Force = 40;
% params.Persistence = 0.4;
% params.ForceSD = 1;
% params.Lc = 28;
% params.LcSD = 3;
% params.numPeaks = 8;
% params.rmsNoise = 5;
% params.nonSpecificForce = 50;
% params.nonSpecificForceSD = 10;
% params.nonSpecificForces = 0;
% params.curveLength = 450;
% params.diffusion = 0.0003;
% [data, Foriginal,rlocs,rpeaks] = generateCurve(params);
% [pks,locs,lcs]=getPeaksLC(data(:,1),data(:,2),1)
% params.Force = 100;
% params.Persistence = 0.4;
% params.ForceSD = 1;
% params.Lc = 28;
% params.LcSD = 3;
% params.numPeaks = 8;
% params.rmsNoise = 5;
% params.nonSpecificForce = 50;
% params.nonSpecificForceSD = 10;
% params.nonSpecificForces = 0;
% params.curveLength = 450;
% params.diffusion = 0.0003;
% [data2, Foriginal,rlocs,rpeaks] = generateCurve(params);
% [pks2,locs2,lcs2]=getPeaksLC(data2(:,1)+80,data2(:,2),1)



numI27 =3;
params.Force = sort([60 60 60 60 200*ones(1,numI27)]);
params.ForceSD = 7.5;
params.Lc = sort([40 60 110 160 220:28:220+28*(numI27-1)]);
params.LcSD = 1;
params.Persistence = 0.4;
params.rmsNoise = 10;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 0;
params.curveLength = 450;
params.diffusion = 0.0006;
[data, Foriginal,rlocs,rpeaks] = generateSpecial(params);
[pks,locs,lcs]=getPeaksLC(data(:,1),data(:,2),1)


numI27 =6;
params.Force = sort([60 60 60 60 200*ones(1,numI27)]);
params.ForceSD = 7.5;
params.Lc = sort([40 70 90 140  200:28:190+28*(numI27-1)]);
params.LcSD = 1;
params.Persistence = 0.4;
params.rmsNoise = 10;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 0;
params.curveLength = 450;
params.diffusion = 0.0006;
[data2, Foriginal,rlocs,rpeaks] = generateSpecial(params);
[pks2,locs2,lcs2]=getPeaksLC(data2(:,1),data2(:,2),1)




bins = 0:4:500
[n1] = hist(lcs,bins)
[n2] = hist(lcs2,bins)

plot(bins,n1,bins,n2)

ma = zeros(length(n1),length(n2));
for i=1:length(n1)
    ma(i,:) = ma(i,:) + n1;
end
for i=1:length(n2)
    ma(:,i) = ma(:,i).*n2';
end
ma = ma/sum(sum(ma));
imagesc(ma)

maxVal = 0;
bestI = [];
ma2 = ma;
for i=-round(length(n2)/2):round(length(n2)/2)
    slice = sum(flipud(diag(diag(ma2,i),i)));
    sortedSlice = sort(slice);
    metric = sum(slice); %sum(sortedSlice(end-50:end));

    subplot(2,1,1)
    imagesc(flipud(ma2+10*diag(diag(ma2,i),i)))
    subplot(2,1,2)
    plot(slice)
    title(sprintf('sortedSliceSum: %2.4f',metric))
    if metric > maxVal
        maxVal = metric;
        bestI = [i];
    elseif metric == maxVal
            bestI = [bestI; i];
    end
%      pause(0.05);
end
bestI
mean(bestI)
maxVal

subplot(3,1,1)
plot(data(:,1),data(:,2),data2(:,1)+mean(diff(bins))*mean(bestI),data2(:,2))
subplot(3,1,2)
imagesc(ma)
subplot(3,1,3)
plot(bins,n1,bins+mean(diff(bins))*mean(bestI),n2)

