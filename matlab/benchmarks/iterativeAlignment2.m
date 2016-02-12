close all;
clear a;
a{1}.L=[10 30 90]+30;
a{1}.adj = 0;
a{2}.L=[10 30 90]-20;
a{2}.adj = 0;
a{3}.L=[10 30 90]+20;
a{3}.adj = 0;
a{4}.L=[10 30 90]-10;
a{4}.adj = 0;
a{5}.L=[10 30 90]+10;
a{5}.adj = 0;
a{6}.L=[10 30 90 ]+32;
a{6}.adj = 0;
a{7}.L=[10 30 90 120 150 180 210]+rand(1,1)*10;
a{7}.adj = 0;
a{8}.L=[10 30 90 120 150 180 210]+rand(1,1)*10;
a{8}.adj = 0;
a{9}.L=[10 30 90 120 150 210]-rand(1,1)*20;
a{9}.adj = 0;
drawPieces(a)

lastMovementMean = 0;
lastMovement = [0 0 0];
for kkk=1:1000
    movement = 0;
    iArray = randperm(length(a));
    jArray = randperm(length(a));
    for ii=1:length(a)
        i=iArray(ii);
        newAdj = [];
        for jj=1:length(a)
            j=jArray(jj);
            if j~=i
                [shiftAmount] = getAlignmentDifference(a{j}.L-a{j}.adj,a{i}.L-a{i}.adj);
                newAdj = [newAdj; shiftAmount];  
            end
        end
        if mean(newAdj) > 0
            a{i}.adj = a{i}.adj + 1;
            movement = movement + mean(newAdj);
        elseif mean(newAdj) < 0
            a{i}.adj = a{i}.adj - 1;
            movement = movement - mean(newAdj);
        end
    end
    drawPieces(a)
    pause(0.1)
    lastMovement = [movement lastMovement];
    lastMovement = lastMovement(1:3);
    if abs(mean(lastMovement) - lastMovementMean) < 1
        break
    end
    lastMovementMean = mean(lastMovement);
    disp(lastMovementMean)
%     if abs(movement) == length(a) || movement == 0
%         break
%     end
end