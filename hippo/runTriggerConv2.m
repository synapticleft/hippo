function [phi whiteningMatrix] = runTriggerConv2(X,pos,thresh,phi)
%% run convolutional sparse coding on TEMPORALLY whitened data, 
%% truncate, then bin and render activations. HAVENT GOTTEN THIS TO WORK

ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
peakToPeak = ceil(1250/dec/8);
%%Processing of position information
bounds = [.1 .9];
pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
vel = angVel(pos);vel = vel(:,1);
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,nanmean(pos));
[~,~,c] = svd(pos(~nanInds,1:2),'econ');pos = (c\pos(:,1:2)')';%pos = a;pos(nanInds) = nan;
pos = pos(:,1);
for i = 1:size(pos,2)   
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
%    pos(:,i) = min(pos(:,i),.9999);
end
pos(nanInds) = 0;
vel = filtLow(vel,1250/32,.5);
vela = vel;
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
vel = vel/max(vel);
vel = resample(vel,ratio,1);
pos = resample(pos,ratio,1);
pos = pos(1:size(X,2),:); 
vel = vel(1:size(X,2));
inds = vel > thresh(1);
%inds(pos < bounds(1) | pos > bounds(2)) = 0;
reg = bwlabel(inds);
h = hist(reg,0:max(reg));
a = accumarray(reg'+1,pos,[],@mean);
f = find(h(2:end) < 1250/dec*thresh(2) | a(2:end)' < bounds(1) | a(2:end)' > bounds(2));
inds(ismember(reg,f)) = 0;
reg = bwlabel(inds);
%% whiten X
%X = bsxfun(@minus,X,mean(X,2));
%S = 64;		% time points in original sources 
J = 64;		% number of basis functions for source generation
R = 20;%20;		% number of time points in basis functions generating sources
N = size(X,1);		% number of sources
[X,whiteningMatrix,dewhiteningMatrix,zerophaseMatrix] = whitenLarge(X,N,50,inds);
reg = reg(inds);
randn('seed',1);
rand('seed',1);

Jrows = 4;

save_every = 200;
display_every = 5;
mintype_inf = 'lbfgsb';
%mintype_lrn = 'minimize';
%lrn_searches = 3;
mintype_lrn = 'gd';
opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 0.01, 'cb', @cb_a);
eta = 0.0001;
eta_up = 1.01;
eta_down = 0.99;
eta_log = [];
target_angle = .01;%0.000001;

paramstr = sprintf('%s_J=%03d_R=%03d_N=%03d_%s', 'hippo', J, R, N, datestr(now,30));
update = 1;
if ~exist('phi','var') || isempty(phi)
phi = randn(N,J,R);
else
    phi1 = phi;
    [~,J,R] = size(phi);
    phi = zeros(N,J,R);
%    for j = 1:J
%        phi(:,j,:) =  whiteningMatrix * squeeze(phi1(:,j,:));
%    end
end
% renormalize
for j = 1:J
    phi(:,j,:) = phi(:,j,:) / sqrt(sum(sum(phi(:,j,:).^2)));
end

num_trials = 200;
% if J == 1
%     lambda = [.2 .25];
% else
%     lambda = [.1 3/50];
% end
lambda = [.05 3/50];
%lambda = .6;
for q = 1:20
    totResp = zeros(J,1);
    sparsenet
%    phi = timeshift_phi(phi,totResp);
%     if lambda < 1.2
%         lambda = lambda + .2;
%     end
    if lambda(1) < .2
        lambda(1) = lambda(1) + .05;
    end
end

function [Y,wh,dw,zp] = whitenLarge(X,N,R,inds)
tic;
X = bsxfun(@minus, X, mean(X(:,inds),2));
f = find(inds);
chunks = 5000;
Xcov = zeros(R);
Xtoep = zeros(R,size(X,1)*chunks);
for i = 1:floor(numel(f)/chunks)
    for j = 1:R
        %Xtoep((j-1)*size(X,1)+(1:size(X,1)),:) = X(:,f((i-1)*chunks+(1:chunks))-R/2+j);
        Xtoep(j,:) = reshape(X(:,f((i-1)*chunks+(1:chunks))-R/2+j),[1 size(X,1)*chunks]);
    end
        Xcov = Xcov + Xtoep*Xtoep';
    toc
end
A = Xcov/(floor(numel(f)/chunks)*chunks);
[V,D] = eig(A);
d = diag(D);
lambda = 0;%exp(1);
dsqrtinv = real(d + lambda).^(-0.5);
wh = diag(dsqrtinv)*V';
dw = V*diag((d+lambda).^.5);%sqrt(D);
zp = V*wh;
Y = zeros(N,numel(f));
    for i = 1:R
        Y = Y + zp(R/2,i)*X(:,f-R/2+i);
    end
toc
clear X;