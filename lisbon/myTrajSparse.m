function [phi] = myTrajSparse(fn)%dewhiteningMatrix,phi) whiteningMatrix
%% run convolutional sparse coding, then bin and render activations.

file = dir('*.mat');
load(file(fn).name,'data');
inds = [10:13 15 16];
choice = [data{2:end,7}];
f = choice ~= 3;
X = [];reg = [];
ord = 20;
XCov = zeros(ord*numel(inds));Xy = zeros(2,numel(inds)*ord);
X = [];
for i = 2:size(data,1)
    if f(i-1)
        trialData = reshape([data{i,inds}],[numel(data{i,6}) numel(inds)])';
        trialData = trialData(:,~isnan(data{i,6}));
        %trialData = diff(trialData,[],2);
        tempXX=  [];
        for j = 1:numel(inds)
            tempXX = [tempXX; toeplitz(trialData(j,:),zeros(1,ord))'];
        end
        %;toeplitz(trialData(2,:),zeros(1,ord+1))'];
        tempXX = tempXX(:,ord:end);
        X = [X tempXX];
%        XCov = XCov + tempXX([2:ord+1 ord+3:end],:)*tempXX([2:ord+1 ord+3:end],:)';
%        Xy = Xy + tempXX([1 ord+2],:)*tempXX([2:ord+1 ord+3:end],:)';
        %X = [X trialData];
        reg = [reg (i-1)*ones(1,size(tempXX,2))];
    end
end
%f = find(f);
%w = Xy/XCov;
%X = w*XFull([2:ord+1 ord+3:end],:) - XFull([1 ord+2],:);%XFull([1 ord+2],:);%
X = bsxfun(@minus,X,mean(X,2));
X = bsxfun(@rdivide,X,std(X,0,2));
%% whiten all signals at once
[X,V] = zca2(X);
X = X(1:ord:end,:);
%% whiten X
%X = bsxfun(@minus,X,mean(X,2));
%S = 64;		% time points in original sources 
J = 32;		% number of basis functions for source generation
R = ord;%20;		% number of time points in basis functions generating sources
% if J >= 1
%     if ~exist('dewhiteningMatrix','var') || isempty(dewhiteningMatrix)
%         numSamples = 50000;
%         indsSub = rand(numel(inds),1) < numSamples/sum(inds);
%         [~,whiteningMatrix,dewhiteningMatrix] = whiten(X(:,indsSub' & inds)');
%     else
%         whiteningMatrix = pinv(dewhiteningMatrix);
%     end
% else
%     numSamples = 50000;
%     indsSub = rand(numel(inds),1) < numSamples/sum(inds);
%     dewhiteningMatrix = diag(std(X(:,indsSub' & inds),0,2));
%     whiteningMatrix = diag(1./diag(dewhiteningMatrix));
% end
% X = whiteningMatrix * X;
whiteningMatrix = eye(size(X,1));
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

num_trials = 400;
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
%    phi = timeshift_phi(phi,totResp);
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
