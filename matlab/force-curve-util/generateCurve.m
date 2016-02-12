function [data,Foriginal,rlocs,rpks] = generateCurve(params)
% GENERATE A CURVE
% [data,Foriginal,rlocs,rpeaks] = generateCurve(params)
% data = [x,y] of generatedCurve
% Foriginal is y without noise
%
% params.Force = 100;
% params.ForceSD = 0;
% params.Lc = 28;
% params.LcSD = 0;
% params.Persistence = 0.4;
% params.numPeaks = 12;
% params.rmsNoise = 10;
% params.nonSpecificForce = 50;
% params.nonSpecificForceSD = 10;
% params.nonSpecificForces = 0;
% params.curveLength = 450;
% params.diffusion = 0.0006;
% data = generateCurve(params);


% Force level and SD
Fu = params.Force;
Fu_sd = params.ForceSD;

% Lc
P = params.Persistence;
Lcbaseline=params.Lc;
startLc_sd = params.LcSD;

% Num peaks
num = params.numPeaks;

% RMS Noise
rmsNoise = params.rmsNoise;

% Nonspecific force
Fu_nonspecific = params.nonSpecificForce;
Fu_nonspecific_sd = params.nonSpecificForceSD;



xbaseline=linspace(0,params.curveLength,5000);
[F] = WLC2(P,xbaseline,Lcbaseline*(num+1));
F(F>Fu) = 0;
Fbaseline = F;

%% Nonspecific
if params.nonSpecificForces > 0
    Lc = Lcbaseline;
    lastLc = 0;
    for kkkk=1:3
    [F] = WLC2(0.5+randn(1,1)*.2,xbaseline,40+randn(1,1)*20);
    F(F>Fu_nonspecific + randn(1,1)*Fu_nonspecific_sd) = 0;
    for i=1:length(Fbaseline)
        if F(i)<Fbaseline(i)
            F(i) = 0;
        else
            F(i) = F(i) - Fbaseline(i);
        end
    end
    Fbaseline = Fbaseline + F;
    lastLc = lastLc + Lcbaseline;
    end
end

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



%% Add Brownian noise
% Set configuration parameters
Ddt = params.diffusion;        % Diffusion coefficient times the time step
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
% close all;
% plot(data(:,1),data(:,2))

[pks,locs] = findpeaks(Foriginal,'MinPeakProminence',2,'MinPeakDistance',5,'MinPeakHeight',5);
rpks = pks;
for i=1:length(locs)
    locs(i) = xbaseline(locs(i));
end
rlocs = locs;

end