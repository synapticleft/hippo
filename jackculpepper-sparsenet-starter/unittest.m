%L = 144;
%M = 4*L;
L =  size(Xf,1);
M = 100;
%Mrows = sqrt(L);

%Lsz = sqrt(L);

tol_coef = 0.01;

display_every = 10;
save_every = 0;% 100;
test_every = 100;

%% sparse penalty
lambda = 2;

%% weight decay
gamma = 0.0;

%% kinds of data
datasource = 'lfp';%'images';

%% number of patches in test set
Btest = 100;

mintype_inf = 'batch';%'lbfgsb';

%% init for optimization method
switch mintype_inf
    case 'lbfgsb'
 
        lb  = zeros(1,M); % lower bound
        ub  = zeros(1,M); % upper bound
        nb  = ones(1,M);  % bound type (lower only)
        nb  = zeros(1,M); % bound type (none)
 
        opts = lbfgs_options('iprint', -1, 'maxits', 20, ...
                             'factr', 1e-1, ...
                             'cb', @cb_a);
    case 'batch'
        opts = lbfgs_options('iprint', -1, 'maxits', 20, ...
                             'factr', 1e-1, ...
                             'cb', @cb_a);
end

mintype_lrn = 'gd';

target_angle_init = 0.05;
target_angle = target_angle_init;

paramstr = sprintf('L=%03d_M=%03d_%s',L,M,datestr(now,30));

B = 100;

reinit

eta = 0.1;
num_trials = 2000;

% for B = 10:10:100
%     sparsenet
% end

for target_angle = target_angle_init:-0.01:0.01
    sparsenet
end


