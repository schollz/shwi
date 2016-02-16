
function [bestMatch] = getBestMatch(a,bOriginal)
% returns the number to subtract from b to be closeset to a
% [bestMatch] = getBestMatch(a,b)

% a=[10 30 60 90];
% bOriginal=[31 57 70 88];
posibilities = [];
num = 0;
for ai=1:length(a)
    for bi=1:length(bOriginal)
        bMod =  (bOriginal(bi)-a(ai));
        b = bOriginal - bMod;
        aMatched = getMatches2(reshape(a,1,length(a)),reshape(b,1,length(b)),10);
        if sum(aMatched>0)==0
            continue
        end
        theDiff = 0;
        for i=1:length(a)
            if aMatched(i) > 0
                theDiff = theDiff + sqrt((a(i)-b(aMatched(i)))^2);
            end
        end
        num = num + 1;
        bestMatching{num} = aMatched;
        posibilities = [posibilities; sum(aMatched>0) 1/theDiff bMod num];
    end
end
posibilities = sortrows(posibilities);
bestMatch = bestMatching{posibilities(end,end)};

end

