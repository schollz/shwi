addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
addpath('../peaks/');
addpath('../benchmarks');

params.Force = 40;
params.Persistence = 0.4;
params.ForceSD = 1;
params.Lc = 28;
params.LcSD = 3;
params.numPeaks = 8;
params.rmsNoise = 5;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 0;
params.curveLength = 450;
params.diffusion = 0.0003;
[data, Foriginal,rlocs,rpeaks] = generateCurve(params);
[pks,locs,lcs]=getPeaksLC(data(:,1),data(:,2),1)
params.Force = 100;
params.Persistence = 0.4;
params.ForceSD = 1;
params.Lc = 28;
params.LcSD = 3;
params.numPeaks = 8;
params.rmsNoise = 5;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 0;
params.curveLength = 450;
params.diffusion = 0.0003;
[data2, Foriginal,rlocs,rpeaks] = generateCurve(params);
[pks2,locs2,lcs2]=getPeaksLC(data2(:,1),data2(:,2),1)




numI27 =randi(9,1,1)-1;
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
[pks,locs,lcs]=getPeaksLC(data(:,1),data(:,2),1)


numI27 =randi(9,1,1)-1;
params.Force = sort([60 60 60 200*ones(1,numI27)]);
params.ForceSD = 7.5;
params.Lc = sort([40 90 140 200:28:190+28*(numI27-1)]);
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


    
hist(lcs2,100)
[f2,x2] = ecdf(lcs2)
plot(x2,f2)
hist(lcs,100)
[f1,x1] = ecdf(lcs)
plot(x1,f1,x2,f2)

bestMetric = 1000;
bestShift = 0;
for shift=-100:2:100
    curMetric = 0;
    for i=1:length(x2)
        tmp = abs(x2-(x2(i)+shift));
        [idx idx] = min(tmp);
        val=x2(idx);
        if val > min(x1) && val < max(x1)
            tmp = abs(x1-val);
            [idx idx] = min(tmp);
            metric = abs(f2(i)-f1(idx));
%               disp(sprintf('%2.1f %2.1f %d %d %2.1f',val,x1(idx),shift,i,metric));
            if metric > curMetric
                curMetric = metric;
            end
        end
    end
    if curMetric < bestMetric
        bestMetric = curMetric;
        bestShift = shift;
    end
end
subplot(2,1,1)
plot(x1,f1,x2+bestShift,f2)
subplot(2,1,2)
plot(data(:,1),data(:,2),data2(:,1)+bestShift,data(:,2))

