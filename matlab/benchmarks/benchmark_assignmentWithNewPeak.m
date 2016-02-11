lcSDs = [1 2 3 4 5 10];
deltaLcs = [1 2 4 8 16 32 64];
data = zeros(length(lcSDs),length(deltaLcs));
for ii=1:length(lcSDs)
    for jj=1:length(deltaLcs)

        deltaLc = deltaLcs(jj);
        lcSD = lcSDs(ii);
% Simulate contour length increments
clear r
trueCidx = [];
for i=1:20
    r{i}.L = [ 30 + lcSD.*randn(1,1) 30 + deltaLcs + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)];
    trueCidx(i) = 1;
end
for i=20:40
    r{i}.L = [ 30 + lcSD.*randn(1,1) 100 + lcSD.*randn(1,1)];
    trueCidx(i) = 2;
end
% for i=60:70
%     r{i}.L = [ 10 + 7.*randn(1,1)  30 + 7.*randn(1,1) 170 + 7.*randn(1,1)]
% end

%% Calculate distance matrix
[r distMatrix] = calculatingMatrix(r);
numClusters = 2;
Z = linkage(distMatrix,'ward','euclidean');
close all; figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');
ylabel('Group Number')
cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;
    

%% Analyze efficiency of assignment (for in 
percentWrong = 0;
for i=1:max(trueCidx)
    subCidx = cidx(find(trueCidx==i));
    percentWrong = percentWrong+ length(find(subCidx ~= mode(subCidx)))/length(cidx);
end
disp(percentWrong)
data(ii,jj) = percentWrong;
    end
end
imagesc(data)
xticks = linspace(1, size(data, 2), size(data, 2));
set(gca,'XTick',xticks,'XTickLabel', deltaLcs)
xticks = linspace(1, size(data, 1), size(data, 1));
set(gca,'YTick',xticks,'YTickLabel', lcSDs)
ylabel('Lc SD [nm]')
xlabel('Delta Lc shift [nm]')

