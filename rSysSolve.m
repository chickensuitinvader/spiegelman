function [ times,hists ] = rSysSolve(d,b,gamma,params,print)

% default print behaviour
if ~exist('print','var'); print = false; end;

% system type
mode = params.mode;

% length parameter and initial number
l = params.l;
n0 = params.n0;

% solver options (time)
tlims = params.tlims;
epochs = params.epochs;

% system parameters
k = params.k; % constant
r = params.r; % number of replicators
reject = params.reject; % rejection boundary
kappa = params.kappa; % fitness function parameter

% output structures
hists = cell(length(kappa),length(gamma));
times = cell(length(kappa),length(gamma));

if print
    tic
end

% for each parameter value of the fitness function/transfer function
for p = 1:length(kappa)
    for q = 1:length(gamma)
        odesys = @(t,n) replicationsystem(t,n,l,d,b,k,r,(l>reject)',mode);
        n = n0; % initial condition
        % data structures for current parameters
        hist = [];
        ts = 0;
        for stop = 1:epochs
            % solve
            [t,n] = ode45(odesys,tlims,n);
            % store epoch data
            hist = [hist;n];
            ts = [ts;ts(end)+t];
            % perform transfer
            n = n(end,:);
            if params.mode >= 0
                n = transfer(l,n,gamma(q),kappa(p));
            end
        end
        ts = ts(2:end); % remove leading zero
        % store data
        hists{p,q} = hist;
        times{p,q} = ts;
    end
end

if print
    toc
end

end

