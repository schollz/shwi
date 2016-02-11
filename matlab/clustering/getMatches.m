function [vecDoneOriginal,vecDoneOriginal2,meanDifference] = getMatches(originalvec1,originalvec2)
% [vecDoneOriginal,vecDoneOriginal2] = getMatches([1 2 5.5 7.7 9.8],[45.6 47.8 49.9])
% originalvec1 = [1 2 5.5 7.7 9.8];
% originalvec2 = [45.6 47.8 49.9];
% 
% originalvec1 = [114.5463
%   131.5626
%   132.4888
%   157.9461
% ];
% 
% originalvec2 = [  272.7692
%   291.5840
%   297.7875
%   345.8923
%   370.1781
%   397.8759
%   425.4866
%   451.2396
%   474.1461
%   500.4904];

    vec1 = sort(originalvec1,'descend');
    vec2 = sort(originalvec2,'descend');
    vec1 = reshape(vec1,1,length(vec1));
    vec2 = reshape(vec2,1,length(vec2));

    while length(vec1) < length(vec2)
        vec1 = [vec1 0];
    end
    while length(vec2) < length(vec1)
        vec2 = [vec2 0];
    end

    m = zeros(length(vec1),length(vec2));
    for n=1:length(vec1)
        for p=1:length(vec2)
            m(n,p) = abs(vec1(n)-vec2(p));
        end
    end
    originalm = m;

    totalError = 0;
    vecDone = zeros(length(vec1),1);
    numPeaks = size(vecDone,1);
    while sum(vecDone)<numPeaks*(numPeaks+1)/2
        mn = min(min(m));
        [rows,cols,vals] = find(m==mn);
        if length(rows)>1
            rows=rows(1);
            cols=cols(1);
        end
        if vecDone(rows)==0 && sum(vecDone==cols)==0
            vecDone(rows)=cols;
            totalError = totalError + mn;
        end
        m(rows,cols) = 500000;
    end

    error = 0;
    vecDoneOriginal = zeros(length(originalvec1),1);
    differences = zeros(length(originalvec1),1);
    for i=1:length(vecDone)
        originalvec1_index = find(originalvec1==vec1(i));
        originalvec2_index = find(originalvec2==vec2(vecDone(i)));
        if length(originalvec1_index) > 0 && length(originalvec2_index) > 0
            vecDoneOriginal(originalvec1_index) = originalvec2_index;
            differences(i) = originalvec1(originalvec1_index)-originalvec2(originalvec2_index);
        end
    end

meanDifference = mean(differences(abs(differences)>0));
originalvec1new = originalvec1-meanDifference;
  vec1 = sort(originalvec1new,'descend');
    vec2 = sort(originalvec2,'descend');
    vec1 = reshape(vec1,1,length(vec1));
    vec2 = reshape(vec2,1,length(vec2));

    while length(vec1) < length(vec2)
        vec1 = [vec1 0];
    end
    while length(vec2) < length(vec1)
        vec2 = [vec2 0];
    end

    m = zeros(length(vec1),length(vec2));
    for n=1:length(vec1)
        for p=1:length(vec2)
            m(n,p) = abs(vec1(n)-vec2(p));
        end
    end
    originalm = m;

    totalError = 0;
    vecDone = zeros(length(vec1),1);
    numPeaks = size(vecDone,1);
    while sum(vecDone)<numPeaks*(numPeaks+1)/2
        mn = min(min(m));
        [rows,cols,vals] = find(m==mn);
        if length(rows)>1
            rows=rows(1);
            cols=cols(1);
        end
        if vecDone(rows)==0 && sum(vecDone==cols)==0
            vecDone(rows)=cols;
            totalError = totalError + mn;
        end
        m(rows,cols) = 500000;
    end

    error = 0;
    vecDoneOriginal = zeros(length(originalvec1),1);
    vecDoneOriginal2 = zeros(length(originalvec2),1);
    differences = zeros(length(originalvec1),1);
    for i=1:length(vecDone)
        originalvec1_index = find(originalvec1new==vec1(i));
        originalvec2_index = find(originalvec2==vec2(vecDone(i)));
        if length(originalvec1_index) > 0 && length(originalvec2_index) > 0
            vecDoneOriginal(originalvec1_index) = originalvec2_index;
            vecDoneOriginal2(originalvec2_index) = originalvec1_index;
            
            differences(i) = originalvec1(originalvec1_index)-originalvec2(originalvec2_index);
        end
    end


end


