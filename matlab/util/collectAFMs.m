function [r] = collectAFMs(folder)


minPeakDistance = 25;
minPeakHeight = 25;
bestLength = 4;

allLs = [];
vecnum = 1;
filelist = getAllFiles(folder);
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
end