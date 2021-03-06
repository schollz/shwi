function [shiftAmount] = getAlignmentDifference(a,bOriginal)
% returns the number to subtract from b to be closeset to a
% [shiftAmount] = getAlignmentDifference(a,b)

% a=[10 30 60 90];
% bOriginal=[31 57 70 88];
posibilities = [];
for ai=1:length(a)
    for bi=1:length(bOriginal)
        bMod =  (bOriginal(bi)-a(ai));
        b = bOriginal - bMod;
        try
            aMatched = getMatches2(a,b,10);
        catch
            aMatched = getMatches2(a',b',10);            
        end
        if sum(aMatched>0)==0
            continue
        end
        theDiff = 0;
        for i=1:length(a)
            if aMatched(i) > 0
                theDiff = theDiff + sqrt((a(i)-b(aMatched(i)))^2);
            end
        end
        posibilities = [posibilities; sum(aMatched>0) 1/theDiff bMod ];
    end
end
posibilities = sortrows(posibilities);
shiftAmount = posibilities(end,3);

end

