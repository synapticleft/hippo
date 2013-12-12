function nickData(X,pos,file,W,f)

inds = [.94 1.036]*100000;
elecs = [1 15 58];
pad = 1000;

[~,pos,~,y] = fixPos(pos);

ratio = round(size(X,2)/numel(pos));

X = X(:,4*(inds(1)-pad)+1:4*(inds(2)+pad));
pos = pos(inds(1)-pad+1:inds(2)+pad);
y = y(inds(1)-pad+1:inds(2)+pad);
pos = resample(pos,ratio,1);
y = round(resample(y,ratio,1));
Xf = morFilter(X,8,1250/8);
[u,~,v] = svds(Xf,1);
W1 = bsxfun(@times,W,exp(1i*angle(mean(pinv(W)))).');

Xf = W1(f,:)*Xf;
Xfd = bsxfun(@times,Xf,exp(1i*-angle(v')));
%[Xz,Z] = zca2(Xf);
%Xzd = Z*Xfd;
%X=  [X;Xf;Xfd;Xz;Xzd;pos'];
X = [abs(Xf);angle(Xf);angle(Xfd);pos'];
X(1:size(Xf,1),:) = bsxfun(@rdivide,X(1:size(Xf,1),:),std(X(1:size(Xf,1),:),0,2));
X = [X;y];
X = real(X);
X = X(:,pad*4+1:end-pad*4);
X(end,:) = X(end,:)-min(X(end,:));

csvwrite(file,X);