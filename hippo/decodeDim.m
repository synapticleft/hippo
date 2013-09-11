function m = decodeDim(Xf,pos,dec)
%%check how OLE decoder works for different number of electrodes / PCs

%inds = 1:64;%33:64;%[1:61 63:64];
nbins = 100;
basis.n = 50;
basis.s = 70;
%yhPos = zeros(basis.n,nbins);
%i = 0;
%figure;

% Positions to use for decoding...
pvec = linspace(0,pi,nbins); % track length is mapped 0-2*pi for von-mises/fourier bases
[~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);
%Xf = Xf(inds,:);
pos = pos(1:size(Xf,2),:);%
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);vel = vel/max(vel);
[~,pos,thresh,~] = fixPos(pos);
pos(pos > 1) = 3-pos(pos > 1);
pos = pos*pi;
[u,s] = eig(Xf(:,thresh)*Xf(:,thresh)');
s = abs(s);
u = bsxfun(@times,u,exp(-1i*angle(mean(u))));
v = (u(:,1)\Xf)';
Xf = bsxfun(@times,Xf,exp(1i*angle(v.')));
if dec > 1
    for j = 1:size(Xf,1)
        Xfd(j,:) = decimate(Xf(j,:),dec);
    end
    Xf = Xfd;clear Xfd;
    vel = decimate(vel,dec);
    pos = angle(-decimate(exp(1i*pos),dec))+pi;
    thresh = logical(round(decimate(double(thresh),dec)));
end
rn = randperm(sum(thresh));
f = find(thresh);
cut = floor(numel(f)/2);
trInds = f(rn(1:cut));
teInds = f(rn(cut+1:end));
%trInds = thresh;%(rn(1:floor(sum(thresh)/2)));
%teInds = thresh;%(rn(ceil(sum(thresh)/2):end));
[~,yRbf] = get1Dbasis('vonmises',basis.n,pos,basis.s);
XfOr = (u*sqrt(s))\Xf;
XfOr = bsxfun(@minus,XfOr,mean(XfOr,2));
for i = 1:size(Xf,1)
    r = 1:size(Xf,1);%randperm(size(XfOr,1));
Xf = [real(XfOr(r(1:i),:));imag(XfOr(r(1:i),:))];
Xf = Xf.';
W = (Xf(trInds,:)'*Xf(trInds,:) + eye(size(Xf,2))/1000)\(Xf(trInds,:)'*yRbf(trInds,:));
yhat = Xf*(W*dbasis');
yhat = bsxfun(@rdivide,yhat,sqrt(mean(abs(yhat(teInds,:))).^2));
%yhat = exp(yhat);
%yhat = bsxfun(@rdivide,yhat,sum(yhat,2));
[~,maxpost]=max(abs(yhat)');
%subplot(2,3,3*(i-1)+3);imagesc((hist3([pos(teInds) maxpost(teInds)'],[30 30])));xlabel('predicted');ylabel('actual');
%for j = 1:size(yhat,2)
%    yhPos(j,:) = accumarray([min(nbins,ceil(pos(teInds)/2/pi*nbins+eps))],yhat(teInds,j),[nbins 1],@mean);
%end
%someHists(2*i-1,:) = hist(pos(trInds)/2/pi,linspace(0,1,nbins));
%someHists(2*i,:) = accumarray([min(nbins,ceil(pos(teInds)/2/pi*nbins+eps))],vel(trInds),[nbins 1],@mean);
%subplot(2,3,3*(i-1)+1);imagesc(abs(yhPos));axis off;title(Beh(f,2));
%subplot(2,3,3*(i-1)+2);imagesc(abs(bsxfun(@rdivide,yhPos,sqrt(sum(abs(yhPos).^2,2)))));
err = circ_dist(pos(teInds),maxpost(teInds)'/size(yhat,2)*2*pi);
m(i) = median(abs(err));
end