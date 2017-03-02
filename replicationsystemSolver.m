%% Replication System Dynamics
%script for solving replicationsystem ODE
clearvars -except d b
close all
%% Parameters and Initial Conditions

% Variables
MEAN = 1001; % mean of intial genome lengths 
SIGMA = 100; % standard deviation of initial genome lengths
SPAN = 2.5; % number of SIGMAs away there will be non-zero number initial genome
SPACING = 0.4; % spacing of initial genome lengths
DIST = 4; % number of SPAN*SIGMAs away there will be lengths of interest
N0 = 1000; % initial number of genomes

% paramter structure
params = struct();
params.mode = 2; % system type

% length parameter and initial number
params.l = MEAN-SPAN*SIGMA*DIST:SPACING*SIGMA:MEAN+SPAN*SIGMA*DIST;
a = (params.l >= MEAN-SPAN*SIGMA) & (params.l <= MEAN+SPAN*SIGMA);
%n0 = N0* ones(length(l),1).*(l==MEAN)'; % point distribution at MEAN
%n0 = N0/(length(l)/DIST)*a;   % uniform distribution at t0
params.n0 = N0*normpdf(params.l,MEAN,SIGMA).*a;         % normal distribution at t0

% epoch time span and epochs
params.tlims = [0, 1e2];
params.epochs = 400;

% probability parameter
% [point addtion, point deletion, block addition, block deletion]
ps = [0.01, 0.01, 5e-6, 5e-6];

% replicator parameters
params.k = 1;
params.r = 10;
params.reject = 45;

% transfer parameters
gamma = 0.05:0.05:0.45;
params.kappa = [0,1,2,3,4,5];

%% Calculation of Mutation Matrices

% Point Mutations
if ~exist('d','var') || (sum(size(d)) ~= length(params.l)*2)
    fprintf('Point Mutation Matrix Time: ')
    tic
    d = calcDeltaM(params.l,ps,params.mode);
    fprintf('%.3f s\n', toc);
end

% Block Mutations
if ~exist('b','var') || (sum(size(b)) ~= length(params.l)*2)
    fprintf('Block Mutation Matrix Time: ')
    tic
    b = calcDeltaB(params.l,ps,params.mode);
    fprintf('%.3f s\n', toc);
end

%% Solving ODE
[ts,hists] = rSysSolve(d,b,gamma,params,1);

%% Plotting
theOtherPlottingFunction(ts,hists,params,gamma);
%edit plotReplication.m