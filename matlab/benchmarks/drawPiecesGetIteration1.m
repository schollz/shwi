function drawPiecesGetIteration1(a)
    % Draw
    maxVal = 0;
    minVal = 10000;
    for i=1:length(a)
        for j=1:length(a{i}.L)
            if a{i}.L(j)+a{i}.Ladj > maxVal
                maxVal = a{i}.L(j)+a{i}.Ladj;
            end
            if a{i}.L(j)+a{i}.Ladj < minVal
                minVal = a{i}.L(j)+a{i}.Ladj;
            end
        end
    end
    minVal = minVal -1;
    imageView = zeros(length(a),round(maxVal-minVal));
    for i=1:length(a)
        for j=1:length(a{i}.L)
            imageView(i,round(a{i}.L(j)+a{i}.Ladj-minVal))=1;
        end
    end
    imagesc(imageView)
end