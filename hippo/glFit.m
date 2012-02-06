function [mod,temp,stats] = glFit(pos,resp)
tic;
lags = 6;nAhead = floor(lags/2);
pos(pos == -1) = nan;
%v = diff(pos);
v = angVel(pos);
v = filtLow(v',1250/32,3)';
resp = bsxfun(@rdivide,resp,std(resp));
pos(1,:) = [];
vp12 = resp(:,1).*conj(resp(:,2));vp12 = filtLow(vp12,1250/32,1);
vp1 = resp(:,1).*conj(resp(:,1));vp1 = filtLow(vp1,1250/32,1);
vp2 = resp(:,2).*conj(resp(:,2));vp2 = filtLow(vp2,1250/32,1);
vp11 = resp(1:end-1,1).*conj(resp(2:end,1));vp11 = filtLow(vp11,1250/32,1);
vp11 = [0;vp11];%vp11 = gsorth(vp11);
vp12 = vp12*exp(1i*pi/4);%./vp1
state = [toeplitz(real(vp11),nan*ones(lags,1)) toeplitz(imag(vp11),nan*ones(lags,1)) ...
    toeplitz(real(vp12),nan*ones(lags,1)) toeplitz(imag(vp12),nan*ones(lags,1))];
state(1:lags-1,:) = [];v(1:lags-1,:) = [];
state = circshift(state,[nAhead 0]);
%state = whiten(state,1000);
%[mod,er,stats] = glmfit(state,max(eps,v(:,1)),'poisson');
out = scale(v(:,1));
[mod,er,stats] = glmfit(state,out,'normal','link','logit');
temp = corr(v(:,1),out-stats.resid);
[er/1000 temp]
figure;plot(out);hold all;plot(out-stats.resid);
%figure;imagesc(log(hist3([out out-stats.resid],[100 100])));
toc

function in = scale(in)
in = in + .2;
in = in/2;
in = max(in,eps);
in = min(in,1-eps);

function X = whiten(X,fudgefactor)
X = bsxfun(@minus, X, mean(X));
X = bsxfun(@rdivide,X,std(X));
A = X'*X/size(X,1);
[V,D] = eig(A);
X = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
