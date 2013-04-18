function [allResp,W,allRun] = serverICA(pos,Xf,v,subset)
accumbins = 100;thresh = .05;
pos = fixPos(pos);
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).')];
Xf = Xf(:,~any(isnan(pos')));
v = v(~any(isnan(pos')));
pos = pos(~any(isnan(pos')),:);
vel = angVel(pos);%vel = filtLow(vel(:,1),1250/32,1);
vel = [0; vel(:,1)];
vel = filtLow(vel,1250/32,.5);
vel = vel/max(vel);
inds = vel > thresh;
b = nan*ones(size(pos,1),1);
bounds = [.1 .9];
Xf = Xf(:,inds);
pos = pos(inds,1);
pos = pos-min(pos)+eps;
pos = pos/max(pos);
b(pos < bounds(1)) = -1;b(pos > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
w = watershed(b==0);
w = w-1; 
pos(mod(w,2) ==1 ,1) = -pos(mod(w,2) ==1 ,1) + 2*max(pos(:));
pos = ceil(pos*accumbins);
Xf = whiten(Xf);
%%%%%%%%%
if ~subset
    %if ~exist('W','var')
    W = zeros(100,size(Xf,1),size(Xf,1));
    %end
    allResp = zeros(size(W,1),size(W,2),accumbins*2);
    for i = 1:size(W,1)
        W(i,:,:) = myACMNsym(Xf);
        act = squeeze(W(i,:,:))'*Xf;
        for j = 1:size(act,1)
            allResp(i,j,:) = accumarray([ones(size(pos)) pos],act(j,:),[1 2*accumbins],@mean);
        end
        i
    end
    save('fullICA.mat','W','allResp');
else
    Xf = Xf(end:-1:1,:);
    allResp = zeros(sum(2:size(Xf,1)),accumbins*2);
    allRun = zeros(1,size(allResp,1));
    counter = 0;
    for i = 1:size(Xf,1)
        W{i} = myACMNsym(Xf(1:i,:));
        act = W{i}'*Xf(1:i,:);
        for j= 1:size(act,1)
            allResp(counter+j,:) = accumarray([ones(size(pos)) pos],act(j,:),[1 2*accumbins],@mean);
        end
        allRun(counter+(1:size(act,1))) = i;
        counter = counter + size(act,1);
    end
    save('subsetICA.mat','allRun','allResp','W');
end

function W = myACMNsym(x)
inVal = 2.5;
 maxcounter = 100;
  [n,m] = size(x);
 pC = (x*transpose(x))/m;
alphas = ones(n,1)*inVal;
    W = complex(randn(n),randn(n));
  Wold = zeros(n);
  k=0;
a2 = .01;
   while (norm(abs(Wold'*W)-eye(n),'fro')>(n*1e-4) && k < maxcounter)%&& k < 15*n)
       minAbsCos = min(abs(diag(W' * Wold)));
    meanAbsCos = mean(abs(diag(W' * Wold)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g %.3g \n',k, norm(abs(Wold'*W)-eye(n),'fro'),1-minAbsCos,1-meanAbsCos);
         k = k+1;
         Wold = W;
         yy = W'*x;
         localAlpha = alphas*.5;
         absy = abs(yy).^2;
         u = absy + a2;
         u1 = bsxfun(@times,bsxfun(@power,u,localAlpha-1),localAlpha);
         %u2 = bsxfun(@times,bsxfun(@power,u,localAlpha-2),localAlpha.*(localAlpha-1));
         u2 = bsxfun(@times,u1./u,localAlpha-1);
         gRad =  x*(u1.*conj(yy)).'/size(x,2); %%double-check
         ggg = mean(u2.*absy + u1,2);
         B = pC*bsxfun(@times,conj(Wold),mean(u2 .* conj(yy).^2,2).');%mean(u2.*(conj(yy).^2),2);%
         W = bsxfun(@times,Wold,ggg.')-gRad + B;%B.'*conj(Wold);
         p = alphas;
         absy = abs(yy);
         u = absy+a2;
         up = bsxfun(@power,u,p);
         u = log(u);
         sigP = mean(bsxfun(@power,absy,p),2).^(1./p);
         gp = -(p.^(-2)).*log(p) + p.^(-2) - psi(1+1./p)./(p.^2) + ...
             mean(bsxfun(@times,up.*bsxfun(@minus,u,1./p+log(sigP)),1./((sigP.^p).*p)),2);
         ggp = 2*p.^(-3).*log(p) - 3./(p.^3) + psi(1,1+1./p).*p.^(-4)+2*psi(1+1./p).*p.^(-3) + ...
             mean(bsxfun(@times,up.*bsxfun(@plus,(u.^2 - bsxfun(@times,u,2./p)-bsxfun(@times,u,2*log(sigP))),...
             2*(p.^-2) + 2*log(sigP)./p + log(sigP).^2),1./(p.*sigP.^p)),2);
         p = p - gp./ggp;
         alphas = min(max(p,.2),4);
    [E,D] = eig(W'*W);
    W = W * E * inv(sqrt(D)) * E';
   end
%A = dewhiteningMatrix*W;
%Z = zerophaseMatrix*A;
%W = W' * whiteningMatrix;