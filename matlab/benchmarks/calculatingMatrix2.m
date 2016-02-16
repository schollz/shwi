function [r,distMatrix] = calculatingMatrix2(r)
    rLength = length(r);
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
                aDlc = [reshape(diff(r{i}.L),length(r{i}.L)-1,1); 0];
                bDlc = [reshape(diff(r{j}.L),length(r{j}.L)-1,1); 0];
                error = 0;
                numMatched = 0;
                for k=1:length(aMatched)
                    if aMatched(k) > 0 && aDlc(k) > 0 && bDlc(aMatched(k)) > 0 && aDlc(k) < 100 && bDlc(aMatched(k)) < 100
%                         disp(sprintf('%2.1f (%2.1f) %2.1f (%2.1f)',r{i}.L(k),aDlc(k),r{j}.L(aMatched(k)),bDlc(aMatched(k))))
                        error = error + sqrt((aDlc(k)-bDlc(aMatched(k)))^2);
                        numMatched = numMatched + 1;
                    end
                end
                error = error/numMatched;
                if isnan(error)
                    error = 0
                end
                distMatrix(i,j) = error;
                distMatrix(j,i) = error;
            end

        end
    end
    textprogressbar('done');
end