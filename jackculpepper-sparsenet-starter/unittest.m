L =  size(Xf,1);
M = 256;

tol_coef = 0.01;

display_every = 10;
save_every = 100;% 100;
test_every = 100;

%% sparse penalty
lambda = [.2 .06];%10
%% weight decay
gamma = 0.0;

%% kinds of data
datasource = 'lfp';%'images';

%% number of patches in test set
Btest = 100;
B = 100;
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
            lb  = zeros(1,B*M); % lower bound
            ub  = zeros(1,B*M); % upper bound
            nb  = ones(1,B*M);  % bound type (lower only)
            nb  = zeros(1,B*M); % bound type (none)
        opts = lbfgs_options('iprint', -1, 'maxits', 20, ...
                             'factr', 1e-1, ...
                             'cb', @cb_a);
end

mintype_lrn = 'gd';

target_angle_init = 0.05;
%target_angle = target_angle_init;

paramstr = sprintf('L=%03d_M=%03d_%s',L,M,datestr(now,30));


reinit

eta = 0.5;
num_trials = 400;%400;

% for B = 10:10:100
%     sparsenet
% end

for target_angle = target_angle_init:-0.01:0.01 %
    %lambda = min(1,lambda + .1);
    sparsenet
end


