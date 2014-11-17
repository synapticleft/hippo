function m = decodeDemo(Xf,pos)
%% OLE decoder of position from LFP
%% models the track as a set of tiled radial basis functions
%Xf = hilbert-transformed multichannel LFP data (#channels x #timesteps)
%pos = rat's current position on the linear track (#timesteps x 4)

nbins = 100; %discretization of track
basis.n = 50; %number of radial basis functions
basis.s = 70; %parameter controlling the width of each basis
lambda = 1/1000; %ridge parameter for regularized regression

% Positions to use for decoding...
pvec = linspace(0,pi,nbins); % track length is mapped 0-2*pi for von-mises/fourier bases
[~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);
pos = pos(1:size(Xf,2),:); %make sure length(pos) == length(Xf)
[~,pos,thresh,~] = fixPos(pos); %thresh marks times when rat is moving
pos(pos > 1) = 3-pos(pos > 1);
pos = pos*pi; %model the bidirectional traversal of the track as walking around a ring
[u,s] = eig(Xf(:,thresh)*Xf(:,thresh)');
s = abs(s);
u = bsxfun(@times,u,exp(-1i*angle(mean(u))));
v = (u(:,1)\Xf)';
Xf = bsxfun(@times,Xf,exp(1i*angle(v.'))); %DEMODULATION OF LFP
rn = randperm(sum(thresh));
f = find(thresh);
cut = floor(numel(f)/2);
trInds = f(rn(1:cut)); %training set
teInds = f(rn(cut+1:end)); %testing set
[~,yRbf] = get1Dbasis('vonmises',basis.n,pos,basis.s);
XfOr = (u*sqrt(s))\Xf; %PCA of Xf
XfOr = bsxfun(@minus,XfOr,mean(XfOr,2)); %zero-mean
Xf = [real(XfOr);imag(XfOr)];
Xf = Xf.';
W = (Xf(trInds,:)'*Xf(trInds,:) + eye(size(Xf,2))*lambda)\(Xf(trInds,:)'*yRbf(trInds,:));
yhat = Xf*(W*dbasis');
yhat = bsxfun(@rdivide,yhat,sqrt(mean(abs(yhat(teInds,:))).^2));
[~,maxpost]= max(abs(yhat)'); %position that is most likely occupied according to decoder
err = circ_dist(pos(teInds),maxpost(teInds)'/size(yhat,2)*2*pi); %error
m = median(abs(err)); %median error across time