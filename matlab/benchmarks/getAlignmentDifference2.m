function [shiftAmount] = getAlignmentDifference2(a,bOriginal,aMatched)
% returns the number to subtract from b to be closeset to a
% [shiftAmount] = getAlignmentDifference(a,b)
% a=[10 30 60 90]-1;
% bOriginal=[31 57 70 88];
shiftAmount = 0;
for i=1:length(aMatched)
    if aMatched(i) > 0
        shiftAmount = shiftAmount + (bOriginal(aMatched(i))-a(i));
    end
end
shiftAmount = shiftAmount / sum(aMatched>0);
end

