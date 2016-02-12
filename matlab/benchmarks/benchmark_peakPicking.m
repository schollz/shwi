forceLevel = [];
theForces = [10 12 14 18 24 30 36 40 80];


for theForceI=1:length(theForces)
  
total2 = 0;
tps2 = 0;
fps2 = 0;
fns2 = 0;
total1 = 0;
tps1 = 0;
fps1 = 0;
fns1 = 0;  

for theIterations=1:10
    theForce = theForces(theForceI);

    params.Force = theForce;
    params.ForceSD = 0;
    params.Lc = 28;
    params.LcSD = 0;
    params.Persistence = 0.4;
    params.numPeaks = 12;
    params.rmsNoise = 10;
    params.nonSpecificForce = 50;
    params.nonSpecificForceSD = 10;
    params.nonSpecificForces = 0;
    params.curveLength = 450;
    params.diffusion = 0.0006;
    [data,Foriginal,rlocs,rpks] = generateCurve(params);

    [pks_g,locs_g] = getPeaksSEGM(data(:,1),data(:,2),1);
    klocs = locs_g;
    kpks = pks_g;

    realPeaks = rlocs;
    guessPeaks = locs_g';
    [vecDoneOriginal] = getMatches2(realPeaks,guessPeaks,5);
    truePositives = sum(vecDoneOriginal>0);
    falseNegatives = sum(vecDoneOriginal==0);
    [vecDoneOriginal] = getMatches2(guessPeaks,realPeaks,5);
    falsePositives = sum(vecDoneOriginal==0);
    total1 = total1 + length(realPeaks);
    tps1 = tps1 + truePositives;
    fps1 = fps1 + falsePositives;
    fns1 = fns1 + falseNegatives;



    minPeakDistance = 12;
    minPeakHeight = 11;
    bestLength = 5;
    [peaks] = getPeaksWFFT(data(:,1),data(:,2),minPeakDistance,minPeakHeight,bestLength);
    wpks = peaks(:,2);
    wlocs = peaks(:,1);
    locs_g = peaks(:,1);
    realPeaks = rlocs;
    guessPeaks = locs_g';
    [vecDoneOriginal] = getMatches2(realPeaks,guessPeaks,5);
    truePositives = sum(vecDoneOriginal>0);
    falseNegatives = sum(vecDoneOriginal==0);
    [vecDoneOriginal] = getMatches2(guessPeaks,realPeaks,5);
    falsePositives = sum(vecDoneOriginal==0);
    total2 = total2 + length(realPeaks);
    tps2 = tps2 + truePositives;
    fps2 = fps2 + falsePositives;
    fns2 = fns2 + falseNegatives;

end


disp(sprintf('PARAMETERS'))
% disp(sprintf('Persistence:\t%2.1f nm \nUnfold Force:\t%2.0f pN\nDelta Lc: \t%2.0f nm\nRMS Noise: \t%2.0f pN',P,Fu,Lcbaseline,rmsNoise))
disp(sprintf('\tTPR\tFPR'))
disp(sprintf('WFFT\t%2.0f\t%2.0f',100*tps2/total2,100*fps2/total2))
disp(sprintf('SEGM\t%2.0f\t%2.0f',100*tps1/total1,100*fps1/total1))

forceLevel = [forceLevel; theForce 100*tps2/total2 100*tps1/total1  100*fps2/total2 100*fps1/total1];

end


subplot(2,1,1)
plot(data(:,1),data(:,2),rlocs,rpks,'+',wlocs,wpks,'*',klocs,kpks,'o','MarkerSize',12)
% text(-40,200,sprintf('Persistence:\t%2.1f nm \nUnfold Force:\t%2.0f pN\nDelta Lc: \t%2.0f nm\nRMS Noise: \t%2.0f pN\nDiffusion: %2.4f',P,Fu,Lcbaseline,rmsNoise,Ddt))
legend('Generated','REAL','WFFT','SEGM','location','NorthEastOutside')
axis([-50 450 -200 400])
subplot(2,1,2)
plot(forceLevel(:,1),forceLevel(:,2),'o-',forceLevel(:,1),forceLevel(:,4),'o-',forceLevel(:,1),forceLevel(:,3),'o-',forceLevel(:,1),forceLevel(:,5),'o-')
legend('WFFT TPR','WFFT FPR','SEGM TPR','SEGM FPR','location','NorthEastOutside')
axis([0 200 -5 105])