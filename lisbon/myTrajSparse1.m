function [phi] = myTrajSparse1(fn)%dewhiteningMatrix,phi) whiteningMatrix
%% run convolutional sparse coding, then bin and render activations.
%% uses cwt to orthogonalize temporal dependencies. also tried toeplitz

%file = dir('*.mat');
load(fn,'data');
inds = [39 31];% 31 [15 16];29 30
indStim = 40;
ord = 4;
choice = [data{2:end,7}];
f = choice ~= 3;% & (1:1339 > 315);
for i = 2:size(data,1)
    if ~numel(data{i,inds(1)})
        f(i-1) = 0;
    end
end
X = [];reg = [];
scales = 2.^(1:3);%(1:.5:2.5);%3.5
wname = 'haar';%'morl';%rbio3.1 and gaus1 are alternatives
for i = 2:size(data,1)
    if f(i-1)
        trialData = reshape([data{i,inds}],[numel(data{i,indStim}) numel(inds)])';
        %fs = [find(~isnan(data{i,indStim}),1) find(~isnan(data{i,indStim}),1,'last')];
        trajcwt = [];
        for j =  1:size(trialData,1)
        %    trajcwt = [trajcwt; toeplitz(zeros(1,ord),trialData(j,:))];%[trajcwt; cwt(trialData(j,:),scales,wname)];
        trajcwt = [trajcwt; diff(trialData(j,:))];
        end
        
        trajcwt = trajcwt(:,~isnan(data{i,indStim}));
        reg = [reg (i-1)*ones(1,size(trajcwt,2))];
        X = [X trajcwt];
        %stimDat = data{i,6} - center_ILD;
        %stimcwt = cwt(stimDat(~isnan(data{i,6})),scales,wname);
        %Xor = [Xor traj(:,~isnan(data{i,6}))];%(f)];
        %stim = [stim stimcwt];
        %stimOr = [stimOr stimDat(~isnan(data{i,6}))];
        %trialNum(i) = size(X,2);
    end
end
f = find(f);
X = zscore(X,0,2);
%X = whiten(X);
%% whiten X
%S = 64;		% time points in original sources 
J = numel(inds);%numel(inds);%1;%		% number of basis functions for source generation
R = 12;		% number of time points in basis functions generating sources

%whiteningMatrix = eye(size(X,1));
N = size(X,1);		% number of sources
randn('seed',1);
rand('seed',1);

Jrows = 4;

save_every = 200;
display_every = 5;
mintype_inf = 'lbfgsb';
%mintype_lrn = 'minimize';
%lrn_searches = 3;
mintype_lrn = 'gd';
%opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 0.01, 'cb', @cb_a);
%       opts.m      Number of limited-memory vectors to use in the algorithm
%                       Try 3 <= m <= 20. (default: 5 )
%opts.factr  Tolerance setting (see this source code for more info)
%                       (default: 1e7 ). This is later multiplied by machine epsilon
%       opts.pgtol  Another tolerance setting, relating to norm(gradient,Inf)
%                       (default: 1e-5)
%       opts.maxTotalIts    How many iterations to allow, including linesearch iterations
%                       (default: 5000)
%       opts.printEvery     How often to display information (default: 1)
%       opts.errFcn 
%opts.maxIts = 30;% (default: 100)

eta = 0.0001;
eta_up = 1.01;
eta_down = 0.99;
eta_log = [];
snrs = [];
energies = [];
target_angle = .01;%0.000001;

paramstr = sprintf('%s_J=%03d_R=%03d_N=%03d_%s', 'hippo', J, R, N, datestr(now,30));
update = 1;
if ~exist('phi','var') || isempty(phi)
phi = randn(N,J,R);
else
    phi1 = phi;
    [~,J,R] = size(phi);
    phi = zeros(N,J,R);
    for j = 1:J
        phi(:,j,:) =  whiteningMatrix * squeeze(phi1(:,j,:));
    end
end
% renormalize
for j = 1:J
    phi(:,j,:) = phi(:,j,:) / sqrt(sum(sum(phi(:,j,:).^2)));
end

num_trials = 100;
% if J == 1
%     lambda = [.2 .25];
% else
     %lambda = [.1 3/50]; tried l(1) = .1, .05
     %lambda = [.1 3/50]; [.05 3/50];
     lambda = [.5 .15];%%[.1 .5];
% end
%lambda = .5;
%lambda = 3;
for q = 1:20
    totResp = zeros(J,1);
    sparsenet
%    lambda = lambda + .05;
    phi = timeshift_phi(phi,totResp);
%    if lambda < 5
%        lambda = lambda + .5;
%    else
if q > 10
        target_angle = target_angle * .9;
    end
    %if lambda < 3
       %lambda = lambda + .5;
    %end
%     if J == 1
%         if lambda(1) < .8
%             lambda(1) = lambda(1) + .2;
%         end
%     else
%     if lambda(1) < .5
%         lambda(1) = lambda(1) + .1;
% %     else
% %         target_angle = target_angle * 0.9;
%     end
%     end
end
