function drawPieces(a)
    % Draw
    maxVal = 0;
    minVal = 10000;
    for i=1:length(a)
        for j=1:length(a{i}.L)
            if a{i}.L(j)-a{i}.adj > maxVal
                maxVal = a{i}.L(j)-a{i}.adj;
            end
            if a{i}.L(j)-a{i}.adj < minVal
                minVal = a{i}.L(j)-a{i}.adj;
            end
        end
    end
    minVal = minVal -1;
    imageView = zeros(length(a),round(maxVal-minVal));
    for i=1:length(a)
        for j=1:length(a{i}.L)
            imageView(i,round(a{i}.L(j)-a{i}.adj-minVal))=1;
        end
    end
    imagesc(imageView)

end