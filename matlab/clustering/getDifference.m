function [meanDifference] = getDifference(originalvec1,originalvec2,vecDoneOriginal1)
% [meanDifference] = getDifference(originalvec1,originalvec2,vecDoneOriginal1)
%   Returns the meanDifference, the value to translate originalvec2 to
%   originalvec1 position
% 
% originalvec1 and originalvec2 come from output of getMatches

    differences = zeros(length(vecDoneOriginal1),1);
    for i=1:length(vecDoneOriginal1)
        if (vecDoneOriginal1(i)>0)
            differences(i) = originalvec1(i) - originalvec2(vecDoneOriginal1(i));
        end
    end
    meanDifference = mean(differences(abs(differences)>0));
end
