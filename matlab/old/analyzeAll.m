clear all;
close all;
warning('off','all')
warning

minPeakDistance = 25;
minPeakHeight = 25;
bestLength = 4;


figure(12)
allLs = [];
vecnum = 1;
filelist = getAllFiles('./');
j = 0;
textprogressbar('collecting files: ');
for i=1:length(filelist)
    textprogressbar(i/length(filelist)*100)
    if  length(findstr('.afm',char(filelist(i))))>0
       data=importdata(char(filelist(i)));
       peaks = analyzeCurve(data(:,1),data(:,2),minPeakDistance,minPeakHeight,bestLength);
       if length(peaks)>0
           j = j + 1;
           rec.file=sprintf('%s',char(filelist(i)));
           rec.x = data(:,1);
           rec.y = data(:,2);
           rec.xPeaks = peaks(:,1);
           rec.yPeaks =  peaks(:,2);
           try
               rec.L=zeros(size(rec.xPeaks));
               for k=1:length(rec.xPeaks)
                    x=rec.xPeaks(k);
                    F=rec.yPeaks(k);
                    kT=4.1;
                    P=0.4;
                    myfun = @(L) F*P/kT-1/(4*(1-x/L).^2)+1/4-x/L;
                    L0=fzero(myfun,round(x)*1.5);
                    rec.L(k) = L0;
               end
               peaks = 1:min(find(rec.yPeaks>150));
               peaks = 1:length(rec.L);
               initL = getInitialContourLength(rec.x,rec.y);
        %            if initL<0
        %                x=rec{k}.x;
        %                y=rec.y;
        %                plot(x,y);
        %                pause(1);
        %            end
               Lcs = [initL; unique(rec.L(peaks))];
               if length(peaks)>1 && length(find(Lcs==Inf | isnan(Lcs)))==0 && initL > 0
                    r{vecnum}.L  = Lcs;
                    r{vecnum}.F = rec.yPeaks(peaks);
                    r{vecnum}.x = rec.x;
                    r{vecnum}.y = rec.y;
                    r{vecnum}.xPeaks = rec.xPeaks(peaks);
%                     close all;
%                     plot(r{vecnum}.x,r{vecnum}.y,'-',r{vecnum}.xPeaks,r{vecnum}.F,'kv','markerfacecolor',[1 0 0],'markeredgecolor',[1 0 0])
%                     pause(0.5)
                    vecnum = vecnum + 1;
               end
           catch
           end
       end
    end
end
textprogressbar('done');


% %% Analyze initial length2
% l1s = [];
% peaks = [];
% allInitialLengths = [];
% allInitialLengthsError = [];
% allNumPeaks =[];
% for i=1:length(r)
%     peaknum =  length(r{i}.F)-1;
%     l1 = r{i}.L(2)-r{i}.L(1);
%     if (l1>0 & l1<54)   
%     peaks = [peaks;  peaknum];
%     l1s = [l1s; l1];
%     end
% end
% plot(peaks,l1s,'o')
% for peaknum=3:8
%     for foocount = 1:100
%         [mu_e,sigma_e,p_e]=gaussian_mixture_model(l1s(peaks==peaknum)',2);
%         [m,mind] = max(p_e);
%         mu_best = mu_e(mind);
%         sigma_best =  sigma_e(mind);
%         if isnan(mu_best) && isnan(sigma_best)
%             [mu_e,sigma_e,p_e]=gaussian_mixture_model(l1s(peaks==peaknum)',3);
%             [m,mind] = max(p_e);
%             mu_best = mu_e(mind);
%             sigma_best =  sigma_e(mind);
%         end
%         if ~isnan(mu_best) && ~isnan(sigma_best)
%             goodind = find(p_e>0.1);
%             allInitialLengths = [allInitialLengths; mu_e(goodind);];
%             allInitialLengthsError = [allInitialLengthsError; sigma_e(goodind);];
%             allNumPeaks = [allNumPeaks; ones(length(goodind),1)*peaknum];
%             break
%         end
%     end
% end
% figure(13)
% hold on;
% plot(allNumPeaks,allInitialLengths,'or')


rLength = length(r);
distMatrix = zeros(rLength);
num = 0;
textprogressbar('calculating matrix: ');
for i=1:rLength-1
    for j=i+1:rLength
        if sum(r{j}.L) ~= sum(r{i}.L)
            num = num + 1;
            textprogressbar((num/(rLength^2/2)*100))
            [vecDoneOriginal,vecDoneOriginal2,diff] = getMatches(r{i}.L,r{j}.L);
            r{i}.matches{j} = vecDoneOriginal; % how i matches to j
            r{i}.diff{j} = diff; % vec j  - vec i
            r{j}.matches{i} = vecDoneOriginal2; % how i matches to j
            r{j}.diff{i} = -1*diff; % vec j  - vec i


            maxLength = length(r{i}.L);
            vec1 = r{i}.L;
            vec2 = r{j}.L;
            vecDone = vecDoneOriginal;
            diff = r{i}.diff{j};
            if length(r{j}.L) > length(r{i}.L)
                maxLength = length(r{j}.L);
                vec1 = r{j}.L;
                vec2 = r{i}.L;
                vecDone = vecDoneOriginal2;
                diff = r{j}.diff{i};
            end

            vec2 = vec2 + diff;
            error = 0;
            for ii=1:length(vecDone)
                if (vecDone(ii)==0)
                    error = error + vec1(ii);
                else
                    error = error + abs(vec1(ii) - vec2(vecDone(ii)));
                end
            end

            distMatrix(i,j) = error;
            distMatrix(j,i) = error;
        end

    end
end
textprogressbar('done');

numClusters = 10;

Z = linkage(distMatrix,'ward','euclidean');
close all;
figure;
subplot(1,2,2)
[H,T,perm]=dendrogram(Z,numClusters,'Orientation','right');


cidx = cluster(Z,'MaxClust',numClusters);
cidx = T;

subplot(1,2,1)
r = iterativeAlignment(r,cidx,perm,100);


%% Plot all
for ci=1:max(cidx)
    cind = find(cidx==ci);
    for jj=1:length(cind)
        i = cind(jj);
        figure(10+ci)
        plot(r{i}.x+r{i}.Ladj,r{i}.y)
        hold on;
    end
end


%% Plot one
ci=4;
cind = find(cidx==ci);
for jj=1:length(cind)
    i = cind(jj);
    figure(10+ci)
    plot(r{i}.x+r{i}.Ladj,r{i}.y)
    hold on;
end





%% Analyze initial length
initialLength = zeros(max(cidx),1);
initialLengthsError = zeros(max(cidx),1);
allInitialLengths = [];
allInitialLengthsError = [];
allNumPeaks = [];
numPeak = zeros(max(cidx),1);
for ci=1:max(cidx)
    cind = find(cidx==ci);
    initialLengths = zeros(length(cind),1);
    numPeaks = zeros(length(cind),1);
    for jj=1:length(cind)
        i = cind(jj);
        initL = r{i}.L(2)-r{i}.L(1);
        if isnan(initL)==0 && initL>0 && initL<200
            initialLengths(jj) = initL;
            numPeaks(jj) = length(r{i}.L)-2;
        else
            initialLengths(jj) = 0;
            numPeaks(jj) = 0;
        end
    end
    initLs = initialLengths(find(initialLengths>0 & initialLengths<800));
    initLValues{ci} = initLs;
    initialLength(ci) = mean(initLs);
    initialLengthsError(ci) = 2*std(initLs)/sqrt(length(initLs));
    if length(initLs)>3
        for foocount = 1:100
            [mu_e,sigma_e,p_e]=gaussian_mixture_model(initLs',3);
            [m,mind] = max(p_e);
            mu_best = mu_e(mind);
            sigma_best =  sigma_e(mind);
            if isnan(mu_best) && isnan(sigma_best)
                [mu_e,sigma_e,p_e]=gaussian_mixture_model(initLs',2);
                [m,mind] = max(p_e);
                mu_best = mu_e(mind);
                sigma_best =  sigma_e(mind);
            end
            if ~isnan(mu_best) && ~isnan(sigma_best)
                initialLength(ci) = mu_best;
                initialLengthsError(ci) = sigma_best;
                goodind = find(p_e>0.2);
                allInitialLengths = [allInitialLengths; mu_e(goodind);];
                allInitialLengthsError = [allInitialLengthsError; sigma_e(goodind);];
                allNumPeaks = [allNumPeaks; ones(length(goodind),1)*mean(numPeaks(numPeaks>0))]
                break
            end
        end
    end
    
    numPeak(ci) = mean(numPeaks(numPeaks>0));
end
figure;
errorbar(numPeak,initialLength,initialLengthsError,'o')

load 8i27test.mat
hold on;
errorbar(numPeak2,initialLength2,initialLengthsError2,'or')


hist([initialLength2./numPeak2; initialLength./numPeak])

vals = initialLength./numPeak;
plot(numPeak,vals,'o')
for i=1:length(vals)
    div = 1
    while vals(i)/div>8
        div = div + 1;
    end
    vals(i) = vals(i) / div;
end
plot(numPeak,vals,'o')

% figure;
% errorbar(allNumPeaks,allInitialLengths,allInitialLengthsError,'o')
% 


% numPeak2 = numPeak;
% initialLength2 = initialLength;
% initialLengthsError2 = initialLengthsError;
% save('8i27test.mat','numPeak2','initialLength2','initialLengthsError2')

numm = 0;
clear masterHist
for ci=1:max(cidx)
    if length(initLValues{ci})>10
        numm = numm + 1
        masterHist{numm} = initLValues{ci}
    end
end
figure;
nhist(masterHist)

nn=1
hist(masterHist{nn})
[mu_e,sigma_e,p_e]=gaussian_mixture_model(masterHist{nn}',3)
[m,mind] = max(p_e);
mu_e(mind)
sigma_e(mind)

