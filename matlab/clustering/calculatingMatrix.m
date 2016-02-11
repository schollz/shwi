function [r,distMatrix] = calculatingMatrix(r)
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
end