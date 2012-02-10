function [mod,temp,stats,outHat] = glFit(pos,resp,inds)
tic;
if ~exist('inds','var')
    inds = 1:size(resp,1);
end
link = 'identity';%'logit';%
lags = 3;nAhead = 0;%floor(lags/2);
pos(pos == -1) = nan;
v = sqrt(sum(diff(pos),2).^2);
%v = diff(pos);
%v = angVel(pos);
v = filtLow(v',1250/32,1)';
v = v(inds,:);resp = resp(inds,:);
%resp = bsxfun(@rdivide,resp,std(resp));
vp12 = resp(:,1).*conj(resp(:,2))./abs(resp(:,1));vp12 = filtLow(vp12,1250/32,1);
vp1 = resp(:,1).*conj(resp(:,1))./abs(resp(:,1));vp1 = filtLow(vp1,1250/32,1);
vp2 = resp(:,2).*conj(resp(:,2));vp2 = filtLow(vp2,1250/32,1);
vp11 = resp(1:end-1,1).*conj(resp(2:end,1))./abs(resp(2:end,1));vp11 = filtLow(vp11,1250/32,1);
vp11 = [0;vp11];vp11 = gsorth(vp11);
%state = [toeplitz(real(vp11),nan*ones(lags,1)) toeplitz(imag(vp11),nan*ones(lags,1)) ...
%   toeplitz(real(vp12),nan*ones(lags,1)) toeplitz(imag(vp12),nan*ones(lags,1))];
state = [toeplitz(real(vp11),nan*ones(lags,1)) toeplitz(imag(vp11),nan*ones(lags,1))];% toeplitz(imag(vp11).*real(vp11),nan*ones(lags,1)) ];%
%state = [toeplitz(real(vp12),nan*ones(lags,1)) toeplitz(imag(vp12),nan*ones(lags,1))];% toeplitz(real(vp12).*imag(vp12),nan*ones(lags,1))]; %toeplitz(vp1,nan*ones(lags,1))];
%state = toeplitz(vp1,nan*ones(lags,1));
state(1:lags-1,:) = [];v(1:lags-1,:) = [];
state = circshift(state,[nAhead 0]);
state = bsxfun(@minus,state,mean(state));
state = bsxfun(@rdivide,state,std(state));
[state dwm] = whiten(state,0);
out = scale(v(:,1));
train = 1:numel(out);test = 1:numel(out);%1:floor(numel(out)/2);test = (1+floor(numel(out)/2)):numel(out); 
[mod,~,stats] = glmfit(state(train,:),out(train),'normal','link',link);

outHat = glmval(mod,state(test,:),link);
temp = corr(max(eps,out(test)),outHat)
mod = [mod(1); dwm*mod(2:end)];
figure;plot(out(test));hold all;plot(outHat);
toc

function a = gsorth(a)
a=complex(real(a),imag(a)-real(a)*diag(sum(real(a).*imag(a))./sum(real(a).^2)));

function in = scale(in)
%in = in + .2;
in = in/10;
in(in<=0) = eps;
in(in >=1) = 1-eps;
%in = max(in,eps);
%in = min(in,1-eps);

function [X DWM] = whiten(X,sub)
A = X'*X/size(X,1);
[V,D] = eig(A);d = diag(D);
if sub
[~, sind] = sort(d,'descend');
d = d(sind);d = d(1:sub);
V = V(:,sind);V = V(:,1:sub);
end
%D1 = sqrt(diag(diag(D)./(diag(D).^2 + fudgefactor^2)));%D1 = sqrt(diag(1./(d+fudgefactor)));
D1 = diag(sqrt(1./d));
WM = D1*V';
X = (WM*X')';
DWM = V*D1^(-1);
%X = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
