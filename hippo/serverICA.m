function [allMax,W] = serverICA(pos,Xf,v,W)
thresh = [.05 1];accumbins = 50;
dec = 32;
%%Processing of position information
bounds = [.2 .9];
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
vel = [0; vel];
pos = bsxfun(@minus,pos,nanmean(pos));
[~,~,c] = svd(pos(~nanInds,1:2),'econ');pos = (c\pos(:,1:2)')';%pos = a;pos(nanInds) = nan;
pos = pos(:,1);
pos = pos - min(pos) + eps;
pos = pos/(max(pos)+eps);
pos(nanInds) = 0;
vel = filtLow(vel,1250/32,.5);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
%vel = vel(1:size(X,2));
inds = vel > thresh(1);
%% which repetition of rat running
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = watershed(b==0);
f = find(runs == 0);
runs(f) = runs(f-1);
%% which contiguous chunk of data
chunk = bwlabel(inds);
h = hist(chunk,0:max(chunk));
a = accumarray([ones(size(chunk)); chunk+1]',pos,[],@mean);
f = find(h(2:end) < 1250/dec*thresh(2) | a(2:end) < bounds(1) | a(2:end) > bounds(2));
inds(ismember(chunk,f)) = 0;
%chunk = bwlabel(inds);
dPos = [0; diff(pos)];
a = accumarray(runs'+1,dPos,[],@mean);
f = find(a(2:end) > 0);
pos(ismember(runs,f)) = 2-pos(ismember(runs,f));pos = pos/2;
posd = floor(pos*accumbins*2)+1;posd = min(2*accumbins,max(1,posd));
posd = posd(inds);v = v(inds,1);
%runs = ceil(runs/2);
%%%%%%%%%
% warning off all;
% dec = 1;
% bounds = [.1 .9];
% pos(pos == -1) = nan;
% reject = 0;
% for i = 1:size(pos,2)
%     reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
% end
% pos(reject,:) = nan;
% if size(Xf,2) < size(pos,1)
%     pos = pos(1:size(Xf,2),:);
% end
% for i = 1:size(pos,2)
%     nanInds = find(~isnan(pos(:,i)));
%     pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
% end
% nanInds = isnan(pos(:,1));
% if size(pos,2) > 2
%     nanInds = nanInds | isnan(pos(:,3));
% end
% pos = pos(~nanInds,:);Xf = Xf(:,~nanInds);%v = v(~nanInds,:);sp = sp(:,~nanInds);
% if dec > 1
%     for i = 1:4
%         posd(:,i) = decimate(pos(:,i),dec);
%     end
%     pos = posd;clear posd;
% end
% pos = bsxfun(@minus,pos,mean(pos));
% % [a,~,~] = svd(pos(:,1:2),'econ');pos = a;
% % for i = 1:2    
% %     pos(:,i) = pos(:,i) - min(pos(:,i));
% %     pos(:,i) = pos(:,i)/(max(pos(:,i)));
% %     pos(:,i) = min(pos(:,i),.9999);
% % end
% %Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).')];
% %Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1)))');
% if dec > 1
% Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
% for i = 2:size(Xf,1)
%     Xfd(i,:) = decimate(Xf(i,:),dec);
% end
% Xf = Xfd;clear Xfd;
% end
% vel = angVel(pos);
% vel = [0; vel(:,1)];
% vel = filtLow(vel,1250/32/dec,1);
% vel = vel/max(vel);
% inds = vel > thresh;
if 0
if ~exist('W','var')
    A = zeros(100,64,63);
    W = zeros(100,63,64);
    Z = zeros(100,64,63);
    al = zeros(100,63);
    for i = 1:100
        [A(i,:,:),W(i,:,:),Z(i,:,:),al(i,:)] = ACMNsym(Xf(:,inds),'mle_circ');
    end
end
allMax = zeros(size(W,1),accumbins*2);
for i = 1:size(W,1)
    act = squeeze(W(i,:,:))*Xf(:,inds);
    act = bsxfun(@times,act,v.');
    acPos = zeros(size(act,1),2*accumbins);
    for j = 1:size(act,1)
        acPos(j,:) = accumarray([ones(size(posd)) posd],act(j,:),[1 2*accumbins],@mean);
    end
    [mxVal,mxInd] = max(abs(acPos)');
    [mxVal,srt] = sort(mxVal,'ascend');
    mxInd = mxInd(srt);
    allMax(i,mxInd) = mxVal;
end
else
    [X,wh] = myWhiten(Xf(:,inds));
    allMax = zeros(size(Xf,1)-1,accumbins*2);
    for i = size(Xf,1)-1:-1:1
        W{i} = myACMNsym(X(i:end,:));
        act = W{i}'*X(i:end,:);
        act = bsxfun(@times,act,v.');
        acPos = zeros(size(act,1),2*accumbins);
        for j= 1:size(act,1)
            acPos(j,:) = accumarray([ones(size(posd)) posd],act(j,:),[1 2*accumbins],@mean);
        end
        [mxVal,mxInd] = max(abs(acPos)');
        [mxVal,srt] = sort(mxVal,'ascend');
        mxInd = mxInd(srt);
        allMax(i,mxInd) = mxVal;
        imagesc(allMax);drawnow;
    end
end

function [x,wh] = myWhiten(Xin)
[Ex, Dx] = eig(cov(Xin'));
% d = flipud(diag(Dx));
% cumVar = sum(d);
% maxLastEig = sum(cumsum(d)/cumVar < .9999999)
% Dx = Dx(end-maxLastEig+1:end,end-maxLastEig+1:end);
% Ex = Ex(:,end-maxLastEig+1:end);
% factors = diag(Dx);
% noise_factors = ones(size(Dx,1),1);
% rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .999999)
% noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
% Dx = diag(factors./noise_factors);
wh = sqrt(inv(Dx)) * Ex';
x = wh * Xin;
% dewhiteningMatrix = Ex * sqrt (Dx);
% %noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
% %D = diag(flipud(1./noise_factors));
% %zerophaseMatrix = E*inv (sqrt (D))*E';
% rolloff_ind = 2;
% noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind)));
% %D = diag(flipud(1./noise_factors));
% zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';%inv (sqrt (D))*E';

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