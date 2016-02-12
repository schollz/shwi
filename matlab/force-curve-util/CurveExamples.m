
% Full Luciferase with Nonspecific
params.Force = [20 30 40 200 200 200 200 200];
params.ForceSD = 7.5;
params.Lc = [30 80 130 190 220 250 280 310];
params.LcSD = 1;
params.Persistence = 0.4;
params.rmsNoise = 10;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 1;
params.curveLength = 450;
params.diffusion = 0.0006;
generateSpecial(params,1)


% Only last two peaks of Lucferase, without Nonspecific
params.Force = [30 40 200 200 200 200 200];
params.ForceSD = 5;
params.Lc = [80 130 190 220 250 280 310];
params.LcSD = 1;
params.Persistence = 0.4;
params.rmsNoise = 10;
params.nonSpecificForce = 50;
params.nonSpecificForceSD = 10;
params.nonSpecificForces = 0;
params.curveLength = 450;
params.diffusion = 0.0006;
generateSpecial(params,1)

