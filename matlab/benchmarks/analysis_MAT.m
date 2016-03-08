close all;
clear all;
addpath('../util');
addpath('../force-curve-util/');
addpath('../clustering/');
addpath('../peaks/');
addpath('../benchmarks');

%% Use a real data
i = 0;
clear filelist
clear r
params.Persistence = 0.4;

fsMain = 'D:\Marszalek Lab\Force Curves\AllYeastPGK\';

fs = {'yPGKCCdumped'}

for jasldf=1:length(fs)
    filelist1 = getAllFiles([fsMain fs{jasldf}]);
    for num=1:length(filelist1)
        if ~isempty(findstr(filelist1{num},'.mat'))==1 
            load(filelist1{num})
            for m=1:length(experiment.recording)
                i = i + 1;
                try
                    r{i}.x = experiment.recording{m}.x;
                    r{i}.y = experiment.recording{m}.y;
                    r{i}.file = experiment.recording{m}.file;
                    r{i}.name = experiment.recording{m}.file;
                    r{i}.dT = 0;

                    r{i}.xPeaks = experiment.recording{m}.xPeaks;
                    r{i}.F = experiment.recording{m}.yPeaks;
                    r{i}.L= experiment.recording{m}.L;
                    if length(r{i}.L) < 1
                        i = i - 1;
                        disp(sprintf('skipping %d',i))  
                    end
                catch
                    i = i - 1;
                end
            end
        end
    end
end
originalR = r;

r=originalR;
[r distMatrix] = calculatingMatrix(r);

figure(10)
numClusters = 7;
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
group = cidx;
subplot(1,2,1)
[r,minR,meanDifferenceSDmean] = iterativeAlignment(r,cidx,perm,10,group);
disp(numClusters)
newCriterion = log(mean(meanDifferenceSDmean)) + log(numClusters);
ylabel('Record Number')
xlabel('Contour length')