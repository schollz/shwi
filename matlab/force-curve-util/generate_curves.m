forceLevel = [];
forces = [10 15 17.5 20 25 30 35 42 50 75];
forces = [10 15 17.5 20 25 30 35 42 50 75 100 150 200 250];
forces = [10 15 18.5 21.5 25 30 35 42 50 80 140 200];
for ijk=1:length(forces)

theForce = forces(ijk);
total2 = 0;
tps2 = 0;
fps2 = 0;
fns2 = 0;
total1 = 0;
tps1 = 0;
fps1 = 0;
fns1 = 0;

for asjldkfj=1:10
% Persistance
P = 0.4;

% Force level and SD
Fu = theForce;
Fu_sd = 0;

% Lc
Lcbaseline=28;
startLc_sd = 0;

% Num peaks
num = 12;

% RMS Noise
rmsNoise = 10;

% Nonspecific force
Fu_nonspecific = 50;
Fu_nonspecific_sd = 10;



xbaseline=linspace(0,450,5000);
[F] = WLC2(P,xbaseline,Lcbaseline*(num+1));
F(F>Fu) = 0;
Fbaseline = F;

%% Nonspecific
% Lc = Lcbaseline;
% lastLc = 0;
% for kkkk=1:3
% [F] = WLC2(0.5+randn(1,1)*.2,xbaseline,40+randn(1,1)*20);
% F(F>Fu_nonspecific + randn(1,1)*Fu_nonspecific_sd) = 0;
% for i=1:length(Fbaseline)
%     if F(i)<Fbaseline(i)
%         F(i) = 0;
%     else
%         F(i) = F(i) - Fbaseline(i);
%     end
% end
% Fbaseline = Fbaseline + F;
% lastLc = lastLc + Lcbaseline;
% plot(xbaseline,Fbaseline)
% end

%% Specific
Lc = Lcbaseline*2 + randn(1,1)*startLc_sd;
lastLc = 0;
realPeaks = zeros(num,2);
for kkkk=1:num
[F] = WLC2(P,xbaseline,Lc);
F(F>Fu + randn(1,1)*Fu_sd) = 0;
for i=1:length(Fbaseline)
    if F(i)<Fbaseline(i)
        F(i) = 0;
    else
        F(i) = F(i) - Fbaseline(i);
    end
end
[y,i] = max(F);
realPeaks(kkkk,1) = xbaseline(i);
realPeaks(kkkk,2) = y+Fbaseline(i);
Fbaseline = Fbaseline + F;
lastLc = lastLc + Lcbaseline;
Lc = Lc + Lcbaseline;
if kkkk==num
    Fbaseline(xbaseline>lastLc)=0;
    break
end

end




%% Add noise
Fbaseline = sgolayfilt(Fbaseline,1,5);
xbaseline = [zeros(1,100)+randn(1,100)*2 xbaseline];
Fbaseline = [linspace(-100,0,100) Fbaseline];
Foriginal = Fbaseline;
noise = wgn(length(Fbaseline),1,0)'*rmsNoise/2;
Fbaseline = Fbaseline + noise;

plot(xbaseline,Foriginal)
axis([-100 500 -150 250])
plot(xbaseline,Fbaseline)
axis([-100 500 -150 250])


%% Add Brownian noise
% Set configuration parameters
Ddt = 0.0006;        % Diffusion coefficient times the time step
sqrtDdt = sqrt(Ddt);
beta = 0;           % Inverse temperature
dwPar = 100;          % A parameter to scale the double well potential

% Set simulation parameters
nSteps = length(Fbaseline);    % Total simulation time (in integration steps)
sampleFreq = 1;     % Sampling frequency
sampleCounter = 1;  % Sampling counter

% Set trajectory to follow particle's position 
xTraj = zeros(1,nSteps/sampleFreq);

% Set initial conditions
x = 0;              % Initial position

for step = 1:nSteps

    % Calculate the new x position after applying a conservative and random force.
    % randn return a random number drawn from a normal distribution with
    % mu=0 and standard deviation=1.
    x = x + 1*beta*Ddt + sqrtDdt*wgn(1,1,0)'*rmsNoise/2;

    % Sample
    if mod(step,sampleFreq) == 0
        xTraj(sampleCounter) = x;
        sampleCounter = sampleCounter + 1;
    end

end
Fbaseline = Fbaseline + xTraj;


%% Finish
data = [xbaseline; Fbaseline]';



[pks,locs] = findpeaks(Foriginal,'MinPeakProminence',2,'MinPeakDistance',5,'MinPeakHeight',5);
rpks = pks;
for i=1:length(locs)
    locs(i) = xbaseline(locs(i));
end
rlocs = locs;


[pks_g,locs_g] = findFEpeaks(data(:,1),data(:,2));
kpks = pks_g;
klocs = locs_g;
a=sort([locs]');
b=sort([locs_g;]);
all = [a zeros(size(a));b  ones(size(b));];
all=sortrows(all);
tps = 0;
total = 0;
fps = 0;
fns = 0;
for i=1:size(all,1)
    theDiff = 1000;
    theI = 0;
    if i<size(all,1)
        if all(i+1,2) ~= all(i,2)
            d = abs(all(i+1,1)-all(i,1));
            if d < theDiff
                theDiff = d;
                theI = i+1;
            end
        end
    end
    if i>1
        if all(i-1,2) ~= all(i,2)
            d = abs(all(i-1,1)-all(i,1));
            if d < theDiff
                theDiff = d;
                theI = i-1;
            end
        end
    end
    if theI > 0
       if all(i,2) == 0
           if theDiff < 5 
               tps = tps + 1;
           else
               fns = fns + 1;
           end
       end
       if all(i,2) == 1
           if theDiff > 5
               fps = fps + 1;
           end
       end
    else
        if all(i,2) == 0
            fns = fns + 1;
        end
        
        if all(i,2) == 1
            fps = fps + 1;
        end
        
    end
   
end
total = total + length(a);
total1 = total1 + total;
tps1 = tps1 + tps;
fps1 = fps1 + fps;
fns1 = fns1 + fns;




minPeakDistance = 12;
minPeakHeight = 11;
bestLength = 5;
[peaks] = analyzeCurve(data(:,1),data(:,2),minPeakDistance,minPeakHeight,bestLength);
wpks = peaks(:,2);
wlocs = peaks(:,1);
locs_g = peaks(:,1);
a=sort([locs]');
b=sort([locs_g;]);
all = [a zeros(size(a));b  ones(size(b));];
all=sortrows(all);
tps = 0;
total = 0;
fps = 0;
fns = 0;
for i=1:size(all,1)
    theDiff = 1000;
    theI = 0;
    if i<size(all,1)
        if all(i+1,2) ~= all(i,2)
            d = abs(all(i+1,1)-all(i,1));
            if d < theDiff
                theDiff = d;
                theI = i+1;
            end
        end
    end
    if i>1
        if all(i-1,2) ~= all(i,2)
            d = abs(all(i-1,1)-all(i,1));
            if d < theDiff
                theDiff = d;
                theI = i-1;
            end
        end
    end
    if theI > 0
       if all(i,2) == 0
           if theDiff < 5 
               tps = tps + 1;
           else
               fns = fns + 1;
           end
       end
       if all(i,2) == 1
           if theDiff > 5
               fps = fps + 1;
           end
       end
    else
        if all(i,2) == 0
            fns = fns + 1;
        end
        
        if all(i,2) == 1
            fps = fps + 1;
        end
        
    end
   
end
total = total + length(a);
total2 = total2 + total;
tps2 = tps2 + tps;
fps2 = fps2 + fps;
fns2 = fns2 + fns;

sum(all)

hold off
plot(data(:,1),data(:,2),rlocs,rpks,'+',wlocs,wpks,'*',klocs,kpks,'o','MarkerSize',12)
legend('Generated','REAL','WFFT','SEGM')
pause(0.1)
end


disp(sprintf('PARAMETERS'))
disp(sprintf('Persistence:\t%2.1f nm \nUnfold Force:\t%2.0f pN\nDelta Lc: \t%2.0f nm\nRMS Noise: \t%2.0f pN',P,Fu,Lcbaseline,rmsNoise))
disp(sprintf('\tTPR\tFPR'))
disp(sprintf('WFFT\t%2.0f\t%2.0f',100*tps2/total2,100*fps2/total2))
disp(sprintf('SEGM\t%2.0f\t%2.0f',100*tps1/total1,100*fps1/total1))

forceLevel = [forceLevel; theForce 100*tps2/total2 100*tps1/total1  100*fps2/total2 100*fps1/total1];

end



subplot(2,1,1)
plot(data(:,1),data(:,2),rlocs,rpks,'+',wlocs,wpks,'*',klocs,kpks,'o','MarkerSize',12)
text(-40,200,sprintf('Persistence:\t%2.1f nm \nUnfold Force:\t%2.0f pN\nDelta Lc: \t%2.0f nm\nRMS Noise: \t%2.0f pN\nDiffusion: %2.4f',P,Fu,Lcbaseline,rmsNoise,Ddt))
legend('Generated','REAL','WFFT','SEGM','location','NorthEastOutside')
axis([-50 450 -200 400])
subplot(2,1,2)
plot(forceLevel(:,1),forceLevel(:,2),'o-',forceLevel(:,1),forceLevel(:,4),'o-',forceLevel(:,1),forceLevel(:,3),'o-',forceLevel(:,1),forceLevel(:,5),'o-')
legend('WFFT TPR','WFFT FPR','SEGM TPR','SEGM FPR','location','NorthEastOutside')
axis([0 200 -5 105])
