theForces = [15 20 25 30 40 60 100];
theLCs = [20 30 40 60];

for method = 1:3
    figure(99)
tpsMatrix = zeros(length(theLCs),length(theForces));
fpsMatrix = zeros(length(theLCs),length(theForces));
for theLcI=1:length(theLCs)

for theForceI=1:length(theForces)
  
    total = 0;
    tps = 0;
    fps = 0;
    fns = 0;

for theIterations=1:50
    theForce = theForces(theForceI);
    theLc = theLCs(theLcI);

    params.Force = theForce;
    params.ForceSD = 0;
    params.Lc = theLc;
    params.LcSD = 0;
    params.Persistence = 0.4;
    params.numPeaks = 8;
    params.rmsNoise = 10;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 0;
    params.curveLength = 450;
    params.diffusion = 0.0006;
    [data,Foriginal,rlocs,rpks] = generateCurve(params);
    
if method == 1
    [pks_g,locs_g] = getPeaksSEGM(data(:,1),data(:,2),1);
    pause(0.2)
    klocs = locs_g;
    kpks = pks_g;
    realPeaks = rlocs;
    guessPeaks = locs_g';
    [vecDoneOriginal] = getMatches2(realPeaks,guessPeaks,10);
    truePositives = sum(vecDoneOriginal>0);
    falseNegatives = sum(vecDoneOriginal==0);
    [vecDoneOriginal] = getMatches2(guessPeaks,realPeaks,10);
    falsePositives = sum(vecDoneOriginal==0);
    total = total + length(realPeaks);
    tps = tps + truePositives;
    fps = fps + falsePositives;
    fns = fns + falseNegatives;
elseif method == 2
    [pks_g,locs_g] = getPeaksLC(data(:,1),data(:,2),1);
    pause(0.2)
    klocs = locs_g';
    kpks = pks_g';
    realPeaks = rlocs;
    guessPeaks = locs_g';
    [vecDoneOriginal] = getMatches2(realPeaks,guessPeaks,10);
    truePositives = sum(vecDoneOriginal>0);
    falseNegatives = sum(vecDoneOriginal==0);
    [vecDoneOriginal] = getMatches2(guessPeaks,realPeaks,10);
    falsePositives = sum(vecDoneOriginal==0);
    total = total + length(realPeaks);
    tps = tps + truePositives;
    fps = fps + falsePositives;
    fns = fns + falseNegatives;
else
    minPeakDistance = 8;
    minPeakHeight = 12;
    bestLength = 5;
    [peaks] = getPeaksWFFT(data(:,1),data(:,2),minPeakDistance,minPeakHeight,bestLength);
    wpks = peaks(:,2);
    wlocs = peaks(:,1);
    locs_g = peaks(:,1);
    realPeaks = rlocs;
    guessPeaks = locs_g';
    [vecDoneOriginal] = getMatches2(realPeaks,guessPeaks,10);
    truePositives = sum(vecDoneOriginal>0);
    falseNegatives = sum(vecDoneOriginal==0);
    [vecDoneOriginal] = getMatches2(guessPeaks,realPeaks,10);
    falsePositives = sum(vecDoneOriginal==0);
    total = total + length(realPeaks);
    tps = tps + truePositives;
    fps = fps + falsePositives;
    fns = fns + falseNegatives;
end

end


tpsMatrix(theLcI,theForceI) = 100*tps/total;
fpsMatrix(theLcI,theForceI) = 100*fps/total;


end

end

figure(method)
subplot(1,2,1)
imagesc(flipud(tpsMatrix))
set(gca,'YTick', 1:length(theLCs),'YTickLabel',fliplr(theLCs));
set(gca,'XTick', 1:length(theForces),'XTickLabel',theForces);
xlabel('Unfolding force [pN]')
ylabel('Contour length increment [nm]')
caxis([0 100])

subplot(1,2,2)
imagesc(flipud(fpsMatrix))
set(gca,'YTick', 1:length(theLCs),'YTickLabel',fliplr(theLCs));
set(gca,'XTick', 1:length(theForces),'XTickLabel',theForces);
xlabel('Unfolding force [pN]')
ylabel('Contour length increment [nm]')
caxis([0 100])

end

