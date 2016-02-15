function [vecDoneOriginal] = getMatches2(originalvec1,originalvec2,cutoff)
% getMatches2 returns the best match within specified distance cutoff
% [vecDoneOriginal] = getMatches2(originalvec1,originalvec2,cutoff)

% originalvec1 = [1.1 10.1 12 3.1];
% originalvec2 = [1 3 10];
% 
% originalvec1 = [69.6739  104.6909  139.6179  174.5449];
% originalvec2 = [23.2246   55.9912   66.6133   79.3959   97.0394  137.0074  137.0074  137.0074  169.5039];
% cutoff = 5;

newVec = [[originalvec1' ones(size(originalvec1))'];[originalvec2' 2*ones(size(originalvec2))']];
newVec = sortrows(newVec);

vecDoneOriginal = zeros(size(originalvec1));
for i=2:size(newVec,1)
    if newVec(i-1,2) ~= newVec(i,2)
        if newVec(i,1)-newVec(i-1,1)<cutoff
           if newVec(i-1,2) == 1
               vecDoneOriginal(find(originalvec1==newVec(i-1,1)))=find(originalvec2==newVec(i,1)); 
           else
               vecDoneOriginal(find(originalvec1==newVec(i,1)))=find(originalvec2==newVec(i-1,1)); 
           end
           
        end
    end
end
end


