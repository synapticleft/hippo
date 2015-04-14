function [phi] = myTrajSparseView(fn,phi)%dewhiteningMatrix,phi) whiteningMatrix
%% run convolutional sparse coding, then bin and render activations.

file = dir('*.mat');
load(file(fn).name,'data');
inds = [15 16];
choice = [data{2:end,7}];
f = choice ~= 3;
X = [];reg = [];
ord = 10;
XCov = zeros(ord*2);Xy = zeros(2,2*ord);
XFull = [];
for i = 2:size(data,1)
    if f(i-1)
        trialData = reshape([data{i,inds}],[numel(data{i,6}) numel(inds)])';
        trialData = trialData(:,~isnan(data{i,6}));
        %trialData = diff(trialData,[],2);
        tempXX = [toeplitz(trialData(1,:),zeros(1,ord+1))';toeplitz(trialData(2,:),zeros(1,ord+1))'];
        tempXX = tempXX(:,ord+1:end);
        XFull = [XFull tempXX];
        XCov = XCov + tempXX([2:ord+1 ord+3:end],:)*tempXX([2:ord+1 ord+3:end],:)';
        Xy = Xy + tempXX([1 ord+2],:)*tempXX([2:ord+1 ord+3:end],:)';
        %X = [X trialData];
        reg = [reg (i-1)*ones(1,size(tempXX,2))];
    end
end
f = find(f);
w = Xy/XCov;
X = w*XFull([2:ord+1 ord+3:end],:) - XFull([1 ord+2],:);%XFull([1 ord+2],:);%
X = bsxfun(@minus,X,mean(X,2));
X = bsxfun(@rdivide,X,std(X,0,2));
J = 32;		% number of basis functions for source generation
R = 20;%20;		% number of time points in basis functions generating sources
whiteningMatrix = eye(size(X,1));
N = size(X,1);		% number of sources
randn('seed',1);
rand('seed',1);

Jrows = 4;
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

paramstr = sprintf('%s_J=%03d_R=%03d_N=%03d_%s', 'hippo', J, R, N, datestr(now,30));
update = 1;
if ~exist('phi','var') || isempty(phi)
phi = randn(N,J,R);
for j = 1:J
    phi(:,j,:) = phi(:,j,:) / sqrt(sum(sum(phi(:,j,:).^2)));
end
% else
%     phi1 = phi;
%     [~,J,R] = size(phi);
%     phi = zeros(N,J,R);
%     for j = 1:J
%         phi(:,j,:) =  whiteningMatrix * squeeze(phi1(:,j,:));
%     end
end
% renormalize


%num_trials = 400;
% if J == 1
%     lambda = [.2 .25];
% else
     %lambda = [.1 3/50]; tried l(1) = .1, .05
     %lambda = [.1 3/50]; [.05 3/50];
     lambda = [.5 .15];%%[.1 .5]
     cols = [1 0 0; 0 1 0];
AFull = [];     
counts = [];
ctrs = [];
useInds = 1:32;
%lambda = 5;
for t = 1:numel(f)
    %% select data
    Xsamp = X(:,reg == f(t));
    %% compute the map estimate
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    a1 = 0;ctr = 1;
%    while ctr <10 %max(abs(a1(:))) < 1 && 
%       subplot(3,3,ctr);
    a0 = randn(J, P);
    opts.x0 = a0(:);
    callF   = @(x) objfun_a_conv(x,Xsamp,phi,lambda);
    a1 = lbfgsb(callF,Inf*ones(1,J*P)',Inf*ones(1,J*P)',opts);
%    imagesc(reshape(a1,J,P));
%    ctr = ctr + 1;
%    end
%    ctrs(t) = ctr;
    a1 = reshape(a1, J, P);a1 = a1(:,R/2+(1:S));
    AFull = [AFull a1];
    figure(1);imagesc(a1);%plot(ctrs);
%    counts(t,1) = max(abs(a1(:)));
%        a0 = randn(J, P);
%    opts.x0 = a0(:);
%    callF   = @(x) objfun_a_conv(x,Xsamp,phi,lambda);
%    a1 = lbfgsb(callF,Inf*ones(1,J*P)',Inf*ones(1,J*P)',opts);
%    counts(t,2) = max(abs(a1(:)));
%	imagesc(a1(:,R/2+(1:S)));
     figure(2);
    for j = 1:32
        subplot(4,8,j);
        inds = reg == f(t);
        temp = XFull([1 ord+2],inds);
        ff = abs(a1(useInds(j),:)) > .5;
        if sum(ff)
        plot(temp(1,min(numel(ff),max(1,find(ff,1)+(-10:10)))),temp(2,min(numel(ff),max(1,find(ff,1)+(-10:10)))),'k--','linewidth',.5);hold all;
        scatter(temp(1,ff),temp(2,ff),[],cols(sign(a1(useInds(j),ff))/2+1.5,:),'filled');hold all;
        end
    end
    drawnow;
%scatter(counts(:,1),counts(:,2));drawnow;
end
