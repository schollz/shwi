
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

fsMain = 'D:\Marszalek Lab\Force Curves\Luc1-450\';

fs = {'160218_PBS2x_MLCT_EINSTEIN',
    '160218_MLCT_PBS2x_DAVINCI'}

for jasldf=1:length(fs)
filelist1 = getAllFiles([fsMain fs{jasldf}]);
for num=1:length(filelist1)
    if ~isempty(findstr(filelist1{num},'.afm'))==1 & sum(findstr(filelist1{num},'_'))>2
        sSplit = strsplit(filelist1{num},'\');
        sSplit = sSplit(end);
        sSplit = strsplit(char(sSplit),'_');
        sSplint = char(sSplit(end));
        dT = str2num(sSplint(1:end-4))/100;
        disp(sprintf('%s %2.2f',char(filelist1{num}),dT))
        if dT > 0 & dT < 20000
            i = i + 1;
            try
                disp(num)
                disp(char(filelist1{num}))
                data = importdata(char(filelist1{num}));
                [origin] = getOrigin(data(:,1),data(:,2),1);
                pause(0.05)
                r{i}.x = data(:,1)-origin(:,1);
                r{i}.y = data(:,2)-origin(:,2);
                clear data
                r{i}.file = sprintf('%d',i);
                r{i}.name = char(filelist1{num});
                r{i}.dT = dT;

                [pks,locs]=getPeaksLC(r{i}.x,r{i}.y,1);
                pause(0.6)
                r{i}.xPeaks = locs;
                r{i}.xPeaks = [r{i}.xPeaks; r{i}.x(end)];
                pks = [pks; r{i}.y(end)];
                r{i}.F = [0];
                r{i}.L=[0];
                peaks = [];
                for k=1:length(r{i}.xPeaks)
                    x=r{i}.xPeaks(k);
                    F=pks(k);
                    try
                        L0 = getLc(params.Persistence,x,F);
                        if ~isnan(L0) && L0 > 0.1
                            r{i}.L = [r{i}.L; L0];
                            r{i}.F = [r{i}.F; pks(k)];
                        end
                    catch
                    end
                end
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
[r distMatrix] = calculatingMatrix2(r);

figure(10)
numClusters = 23;
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
ci=17;
cind = find(cidx==ci);
cm=colormap(jet(length(cind)+1));
num = 0
for jj=1:length(cind)
    num = num + 1;
    i = cind(jj);
    figure(100+ci)
    plot(r{i}.x+r{i}.Ladj,r{i}.y,'color',cm(num,:))
    axis([-350 250 -200 400])
    hold on;
end


