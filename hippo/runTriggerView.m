function posInds = runTriggerView(pos,v,Xf,accumbins,thresh,r,probes,posInds)

bounds = [.1 .9];
pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
%plot(pos(:,1),'r');return;
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
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,mean(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
    posd(:,i) = floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1;
end
offSet = 1;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
     [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xf = [real(Xf);imag(Xf)];
vel = filtLow(vel,1250/32,1);
vel = vel/max(vel);inds = vel > thresh;
[E, D]=pcamat(Xf, 1, numel(posInds), 'off','off');
r1 = pinv(r);
r1 = r1(:,posInds);r = r(posInds,:);
r1 = E*inv(sqrt (D))*E'*r1;
t = r*zscore(Xf,0,2);%(:,inds)
%r1 = r';%pinv(r);
xdim = ceil(sqrt(size(t,1)));ydim = ceil(size(t,1)/xdim);
% %%FOR 1D TRACK
% %%v = filtLow(v.',1250/32,1).';
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
vInterp = zeros(2,size(t,1),max(runs),accumbins(1));
w = watershed(b==0);
w = w-1; %w(w== max(w)) = 0;
%%
for k = 1:2
    runs1 = bwlabel(w>0 & mod(w,2) == k-1 & w <=2*max(runs));%b*((-1)^k)>0);
    inds = runs1 > 0;
    for j = 1:size(t,1)
        vInterp(k,j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) accumbins(1)] ,@mean);
    end
end
t1 = [squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))];
%h1 = figure;
spatial = zeros(size(t1,1),2*accumbins(1));
negs = zeros(size(t1,1),1);
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
posInds = 1:numel(posInds);%find(max(spatial') > 300);%1:size(r,1);%
[~,peakLoc] = max(spatial(posInds,:)');
[~,indLoc] = sort(peakLoc);
posInds = posInds(indLoc);
figure;imagesc(spatial(posInds,:))
h1 = figure;
h2 = figure;
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
sk = ones(1,numel(posInds));
tes = zeros(numel(posInds),max(runs),2*accumbins(1));
ups = zeros(numel(posInds),size(probes,1),size(probes,2)+1);
for i = 1:numel(posInds)
    te = reshape(t1(posInds(i),:),[max(runs) 2*accumbins(1)]);
    if skewness(te(:)) < 0
        te = -te;
        sk(i) = -1;
    end
    figure(h1);subplot(xdim,ydim,i);imagesc(te,[0 max(te(:))]);axis off;%s(temp(indLoc(i)))
    tes(i,:,:) = te;
    u = complex(r1(1:size(Xf,1)/2-1,posInds(i)),r1(size(Xf,1)/2+1:end-1,posInds(i)));%r1(1:size(Xf,1)-1,posInds(i));%
    up1 = probes;
    if exist('probes','var')
    for ii = 1:size(probes,1)
        for j = 1:size(probes,2)
            up1(ii,j) = u(probes(ii,j)+1);%-256
        end
    end
%    up1 = up1(:,[1:4 6 5 8 7]);
    up1 = up1(:,[1:12 14 13 16 15]);
    up1 = [up1(:,1:8) zeros(size(up1,1),1) up1(:,9:16)];
    ups(i,:,:) = up1;
    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(up1,0,1));axis off;
    end
end
%sPlot(bsxfun(@times,t(temp(indLoc),:),sk'));
figure;plot(bsxfun(@times,t(posInds,:),sk')');
%figure;imagesc(complexIm(corr(complex(r(temp(indLoc),1:513),r(temp(indLoc),514:end))'),0,1));
figure;imagesc(complexIm(corr(ups(:,:)'),0,1));
figure;imagesc(squeeze(std(ups)));
superImp(tes);