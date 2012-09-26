function [posInds,t,vel] = runTriggerViewT(pos,v,Xf,accumbins,thresh,r,probes,posInds,r1)

pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);%sp = sp(:,~nanInds);
vel = angVel(pos);
vel = [zeros(1,2); vel(:,1:2)];
for i = 1:2
    vel(:,i) = filtLow(vel(:,i),1250/32,1);
end
veld = vel(:,1:2);
vel = vel(:,1);
vel = vel/max(vel);inds = vel > thresh;
pos = bsxfun(@minus,pos,mean(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
    veld(:,i) = veld(:,i) - min(veld(:,i));
    veld(:,i) = veld(:,i)/max(veld(:,i));
    veld(:,i) = min(veld(:,i),.9999);
    posd(:,i) = floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1;
    veld(:,i) = floor(veld(:,i)*accumbins(min(numel(accumbins),i)))+1;
end
% offSet = 1;
% Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
%     [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
% %  Xf = [bsxfun(@times,Xf,v(:,1).');...
% %    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))].'];
% Xf = [real(Xf);imag(Xf)];
% Xf = zscore(Xf,0,2);%bsxfun(@minus,Xf,mean(Xf,2));
% if ~exist('r1','var')
%     r1 = pinv(r);%r';%
% end
if exist('posInds','var') && ~isempty(posInds)
    r1 = r1(:,posInds);r = r(posInds,:); %% IS THIS RIGHT??
end
% t = r*zscore(Xf,0,2);
t = Xf(posInds,:);%clear Xf;
if 0
    [B M] = size(t);
    opts = lbfgs_options('iprint', -1, 'maxits', 20, ...
        'factr', 1e-1, ...
        'cb', @cb_a);
    %a = phi\X;
    M1 = 10000;
    for i = 1:ceil(size(t,2)/M1)
        ind = (i-1)*M1+(1:M1);ind(ind > size(t,2)) = [];
        M = numel(ind);
        lb  = zeros(1,B*M); % lower bound
        ub  = zeros(1,B*M); % upper bound
        nb  = ones(1,B*M);  % bound type (lower only)
        nb  = zeros(1,B*M); % bound type (none)
        temp = t(:,ind);
        t(:,ind) = reshape(lbfgs(@objfun_a,temp(:),lb,ub,nb,opts,pinv(r),zscore(Xf(:,ind),0,2),2),B,M);
    end
end
% %2d stuff
%
if isempty(posInds)
    posInds = 1:size(r1,2);
end
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
f1 = figure;f2 = figure;
t = bsxfun(@times,t, sign(skewness(t,0,2)));
%t = zscore(t,0,2);
for i = 1:numel(posInds)
%    u = complex(r1(1:size(Xf,1)/2-1,i),r1(size(Xf,1)/2+1:end-1,i));%r1(1:size(Xf,1)-1,posInds(i));%
    u = complex(r1(1:end/2-1,i),r1(end/2+1:end-1,i));
    if exist('probes','var') && ~isempty(probes)
        up1 = probes;
        for ii = 1:size(probes,1)
            for j = 1:size(probes,2)
                up1(ii,j) = u(probes(ii,j)+1);%-256
            end
        end
        %    up1 = up1(:,[1:4 6 5 8 7]);
        up1 = up1(:,[1:12 14 13 16 15]);
        %up1 = diff(up1);
        up1 = [up1(:,1:8) zeros(size(up1,1),1) up1(:,9:16)];
    else
        up1 = reshape(u,[8 8]);
    end
    figure(f1);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
    figure(f2);subplot(xdim,ydim,i);imagesc(imfilter(accumarray(posd,t(i,:),accumbins,@mean,0),fspecial('gaussian',5,1)),[-.1 1]);axis off;
end
sPlot([10*vel';t;abs(v(:,1)')/1000]);
return
% %%FOR 1D TRACK
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
vInterp = zeros(2,size(t,1),max(runs),accumbins(1));
w = watershed(b==0);
w = w-1; %w(w== max(w)) = 0;
%u = zeros(max(w),size(Xf,1));
% for i = 1:max(w)/2
%     [temp,~,~] = svds(Xf(:,(w == 2*i | w == 2*i-1) & inds'),2);
%     u(i,:) = temp(:,2);
% end
%
for k = 1:2
    runs1 = bwlabel(w>0 & mod(w,2) == k-1 & w <=2*max(runs));%b*((-1)^k)>0);
    inds = runs1 > 0;
    for j = 1:size(t,1)
        vInterp(k,j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) accumbins(1)] ,@mean);
    end
end
t1 = [squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))];
%h1 = figure;
spatial = randn(size(t1,1),2*accumbins(1));
for i = 1:size(t1,1)
    temp = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    [~,s,v] = svds(temp,1);
    v = s*v'; 
    if -min(v) > max(v) 
        v = -v; 
    end
    spatial(i,:) = s*v';
%    subplot(xdim,ydim,i);imagesc(temp);axis off;
%    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(reshape(complex(r1(1:32,i),r1(34:65,i)),[8 4]),0,1));axis off;
end
figure;plot(spatial');
if ~exist('posInds','var') || isempty(posInds)
    posInds = find(max(spatial') > 0);
else
   posInds = 1:size(r,1);%
end
spatial = spatial(posInds,:);
[~,peakLoc] = max(spatial');
[~,indLoc] = sort(peakLoc);
%posInds = posInds(indLoc);
spatial = spatial(indLoc,:);
figure;imagesc(spatial);
h1 = figure;
h2 = figure;
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
sk = ones(1,numel(posInds));
tes = zeros(numel(posInds),max(runs),2*accumbins(1));
if exist('probes','var') && ~isempty(probes)
    ups = zeros(numel(posInds),size(probes,1),size(probes,2)+1);
else
    ups = zeros(numel(posInds),8,(size(Xf,1)/2-1)/8);%
end
for i = 1:numel(posInds)
    te = reshape(t1(posInds(i),:),[max(runs) 2*accumbins(1)]);
    if skewness(te(:)) < 0
        te = -te;
        sk(i) = -1;
    end
    figure(h1);subplot(xdim,ydim,i);imagesc(imfilter(te,fspecial('gaussian',5,1)));axis off;%s(temp(indLoc(i))),[0 max(te(:))]
    tes(i,:,:) = te;
    u = complex(r1(1:size(Xf,1)/2-1,posInds(i)),r1(size(Xf,1)/2+1:end-1,posInds(i)));%r1(1:size(Xf,1)-1,posInds(i));%
    if exist('probes','var') && ~isempty(probes)
        up1 = probes;
        for ii = 1:size(probes,1)
            for j = 1:size(probes,2)
                up1(ii,j) = u(probes(ii,j)+1-min(probes(:)));%-256
            end
        end
        %up1 = up1(:,[1:4 6 5 8 7]);
        up1 = up1(:,[1:12 14 13 16 15]);
        up1 = [up1(:,1:8) zeros(size(up1,1),1) up1(:,9:16)];
    else
        up1 = reshape(u,[8 ,(size(Xf,1)/2-1)/8]);
    end
    ups(i,:,:) = up1;
    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
end
t = t(posInds,:);
sPlot([bsxfun(@times,t,sk'); vel']);
superImp(tes,[],1);

%%%%%%%%%%snippet
offSet = 1;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xf = [real(Xf);imag(Xf)];
Xf = zscore(Xf,0,2);
t = squeeze(W(end,:,:))*Xf;
clear Xf;
xc = zeros(4,size(t,1),size(t,1));
for i = 1:size(t,1)
    for j = 1:i-1
        temp = xcorr(t(i,:),t(j,:),1000);
        [a b] = min(temp);
        [c d] = max(temp);
        if abs(a) > abs(c)
            xc(:,i,j) = [a b c d];
        else
            xc(:,i,j) = [-c d -a b];
        end
    end
end