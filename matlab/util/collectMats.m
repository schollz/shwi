function [r] = collectMats(folder,maxPeaksSet)
    addpath('../util')
    allLs = [];
    vecnum = 1;
    fileList12313 = getAllFiles(folder);
    for i=1:length(fileList12313)
        hasMat = findstr(char(fileList12313(i)),'s.mat');
        if length(hasMat)>0
           disp(sprintf('%s',char(fileList12313(i)))) 
           load(char(fileList12313(i)))
           rec = experiment.recording;
           for j=1:length(rec)
               if maxPeaksSet <= 0 || maxPeaksSet>length(rec{j}.xPeaks)
                   maxPeaks = length(rec{j}.xPeaks);
               else
                   maxPeaks = maxPeaksSet;
               end
               rec{j}.L=[];
               peaks = [];
               initL = getOrigin(rec{j}.x,rec{j}.y);
               for k=1:maxPeaks
                    x=rec{j}.xPeaks(k)-initL(1);
                    F=rec{j}.yPeaks(k)-initL(2);
                    kT=4.1;
                    P=0.7;
                    myfun = @(L) F*P/kT-1/(4*(1-x/L).^2)+1/4-x/L;
                    L0=fzero(myfun,round(x)*1.5);
                    if x>0
                    rec{j}.L = [rec{j}.L; L0];
                    peaks = [peaks; k];
                    end
               end
               Lcs = [unique(rec{j}.L)];
               if length(peaks)>1 && length(find(Lcs==Inf | isnan(Lcs)))==0 
                    r{vecnum}.L  = Lcs;
                    r{vecnum}.F = rec{j}.yPeaks(peaks);
                    r{vecnum}.x = rec{j}.x;
                    r{vecnum}.y = rec{j}.y;
                    r{vecnum}.name = char(fileList12313(i));
                    r{vecnum}.file = rec{j}.file;
                    vecnum = vecnum + 1;
               end
           end
        end
    end

end