function [data,Foriginal,rlocs,rpks] = generateSpecial(params,toPlot)
% GENERATE A CURVE
% [data,Foriginal,rlocs,rpeaks] = generateSpecial(params)
% data = [x,y] of generated curve
% Foriginal is y without noise
%
% params.Force = [20 30 40 200 200 200 200 200];
% params.ForceSD = 0;
% params.Lc = [30 80 130 190 220 250 280 310];
% params.LcSD = 0;
% params.Persistence = 0.4;
% params.rmsNoise = 10;
% params.nonSpecificForce = 50;
% params.nonSpecificForceSD = 10;
% params.nonSpecificForces = 0;
% params.curveLength = 450;
% params.diffusion = 0.0006;

if nargin < 2
    toPlot = 0;
end
% Force level and SD
Fu = params.Force;
Fu_sd = params.ForceSD;

% Lc
P = params.Persistence;
Lcbaseline=params.Lc;
startLc_sd = params.LcSD;

% RMS Noise
rmsNoise = params.rmsNoise;

% Nonspecific force
Fu_nonspecific = params.nonSpecificForce;
Fu_nonspecific_sd = params.nonSpecificForceSD;



xbaseline=linspace(0,params.curveLength,5000);
[F] = WLC2(P,xbaseline,params.Lc(end));
F(F>max(Fu)) = 0;
Fbaseline = F;

%% Nonspecific
if params.nonSpecificForces > 0
    Lc = Lcbaseline;
    lastLc = 0;
    for kkkk=1:randi(3,1,1)
    [F] = WLC2(0.5+randn(1,1)*.2,xbaseline,5 + (50-5).*rand(1,1));
    F(F>Fu_nonspecific + 5 + (100-5).*rand(1,1)) = 0;
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
for kkkk=1:length(Fu)-1
    [F] = WLC2(P,xbaseline,Lcbaseline(kkkk)+randn(1,1)*startLc_sd);
    F(F>Fu(kkkk) + randn(1,1)*Fu_sd) = 0;
    for i=1:length(Fbaseline)
        if F(i)<Fbaseline(i)
            F(i) = 0;
        else
            F(i) = F(i) - Fbaseline(i);
        end
    end
    Fbaseline = Fbaseline + F;
%     plot(xbaseline,Fbaseline)
%     pause(1)
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

[pks,locs] = findpeaks(Foriginal,'MinPeakProminence',2,'MinPeakDistance',5,'MinPeakHeight',5);
rpks = pks;
for i=1:length(locs)
    locs(i) = xbaseline(locs(i));
end
rlocs = locs;

if toPlot > 0 
    plot(data(:,1),data(:,2),data(:,1),Foriginal,locs,pks,'o')
end

end