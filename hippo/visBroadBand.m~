function visBroadBand(fname,Xf,v,r,posInds,keepInds)

offSet = 1;
Xf1 = Xf;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
%u1 = mean(Xf,2);u1 = u1(1:end-1);um = conj(mean(u1));u1 = u1*um;
%  Xf = [bsxfun(@times,Xf,v(:,1).');...
%    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))].'];
Xf = [real(Xf);imag(Xf)];

if ~exist('posInds','var') || isempty(posInds)
    posInds = 1:size(r,1);
end
Xf = zscore(Xf,0,2);
t = r(posInds,:)*Xf;
r1 = pinv(r);
r1 = r1(:,posInds);r = r(posInds,:);
%sk = skewness(t,0,2);t = bsxfun(@times,t,sk); r = bsxfun(@times,r,sk);r1 = bsxfun(@times,r1',sk)';

if ~exist('keepInds','var')
    sPlot(t);
    keepInds(1) = input('input start time: ');
    keepInds(2) = input('input end time: ');
end
t = t(:,keepInds(1):keepInds(2));Xf = Xf(:,keepInds(1):keepInds(2));Xf1 = Xf1(:,keepInds(1):keepInds(2));
[B M] = size(t);
opts = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 1e-1,'cb', @cb_a);
lb  = zeros(1,B*M); % lower bound
ub  = zeros(1,B*M); % upper bound
nb  = zeros(1,B*M); % bound type (none)
t1 = reshape(lbfgs(@objfun_a,t(:),lb,ub,nb,opts,r1,Xf,2),B,M);
r1 = complex(r1(1:size(Xf,1)/2-1,:),r1((end-size(Xf,1)/2+1):(end-1),:));

xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
figure;
for i = 1:numel(posInds)
    subplot(xdim,ydim,i);imagesc(complexIm(reshape(r1(:,i),[8 size(r1,1)/8]),0,1));axis off;
end
f1 = figure;subplot(411);sPlot(t1,[],0);
X = getData(fname,1:64,range(keepInds)*32,keepInds(1)*32);
X = zscore(X,0,2);
[u,s,v] = svds(X,1);
whichInd = input('index: ');
col = squeeze(complexIm(r1(:,whichInd),0,1));
figure(f1);
set(gcf,'DefaultAxesColorOrder',col);
subplot(412);plot((1:size(X,2))/32,X');hold all;axis tight;%plot(1:size(X,2),s*v','k','linewidth',2);
Xf = morFilter(X,8,1250);
X = X-u*s*v';
[u,s,v] = svds(Xf,1);u = u*conj(mean(u));
subplot(413);plot((1:size(X,2))/32,real(Xf)');axis tight;
Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');
subplot(414);plot((1:size(X,2))/32,real(bsxfun(@rdivide,Xf,u))');axis tight;
%plot((1:size(X,2))/32,X');hold all;plot((1:size(X,2))/32,s*v','k','linewidth',2);axis tight;
