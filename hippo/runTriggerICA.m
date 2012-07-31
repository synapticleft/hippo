function [r t] = runTriggerICA(pos,v,Xf,accumbins,thresh)
warning off all;
dec = 8;
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);%sp = sp(:,~nanInds);
for i = 1:4
    posd(:,i) = decimate(pos(:,i),dec);
end
pos = posd;clear posd;
vel = angVel(pos);%vel = filtLow(vel(:,1),1250/32,1);
vel = [0; vel(:,1)];
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i))+eps);
    pos(:,i) = floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1;
end
offSet = 1;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
     [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
for i = 2:size(Xf,1)
    Xfd(i,:) = decimate(Xf(i,:),dec);
end
Xf = Xfd;clear Xfd;
Xf = [real(Xf);imag(Xf)];%[abs(Xf); angle(Xf)];
Xf = zscore(Xf,0,2);
vel = filtLow(vel,1250/32/dec,1);
figure;plot(vel);drawnow;
vel = vel/max(vel);inds = vel > thresh;
sum(inds)
t = Xf(:,inds);
[r,~,t] = runica(Xf(:,inds),'pca',50);
xdim = ceil(sqrt(size(t,1)));ydim = ceil(size(t,1)/xdim);
figure;for i = 1:size(t,1)
subplot(xdim,ydim,i);imagesc(accumarray(pos(inds,1:2),t(i,:),[],@mean));
end
% Xf = r*Xf;
%%v = filtLow(v.',1250/32,1).';