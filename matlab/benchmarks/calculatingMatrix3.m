function [r,distMatrix] = calculatingMatrix3(r)
    rLength = length(r);
    for i=1:rLength
        r{i}.L = unique(r{i}.L);
        r{i}.F = unique(r{i}.F);
    end
    textprogressbar('matching peaks: ')
    num = 0;
    for i=1:rLength
        for j=1:rLength
            num = num + 1;
            textprogressbar(num/(length(r)^2)*100)
            if i~=j
                r{i}.bestMatch{j} = getBestMatch(r{i}.L,r{j}.L);
            end
        end
    end
    textprogressbar('done.')
    distMatrix = zeros(rLength);
    num = 0;
    textprogressbar('calculating matrix: ');
    for i=1:rLength-1
        for j=i+1:rLength
            num = num + 1;
            textprogressbar((num/(rLength^2/2)*100))
            if sum(r{j}.L) ~= sum(r{i}.L)
                aMatched = r{i}.bestMatch{j};
                error = 0;
                for ij=1:length(aMatched)
                    if aMatched(ij) == 0
                        error = error + r{i}.consistency(ij);
                    end
                end
                aMatched = r{j}.bestMatch{i};
                for ij=1:length(aMatched)
                    if aMatched(ij) == 0
                        error = error + r{j}.consistency(ij);
                    end
                end
               
                distMatrix(i,j) = error;
                distMatrix(j,i) = error;
            end

        end
    end
    textprogressbar('done');
end