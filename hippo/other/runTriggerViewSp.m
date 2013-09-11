function [posInds,t,tes] = runTriggerViewSp(pos,v,Xf,accumbins,thresh,posInds,shank)
%% some unknown artifact

bounds = [.1 .9];
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
%vel = diff(pos);
vel = [zeros(1,2); vel(:,1:2)];
for i = 1:2
    vel(:,i) = filtLow(vel(:,i),1250/32,1);
end
veld = vel(:,1:2);
vel = vel(:,1);
vel = vel/max(vel);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
inds = vel > thresh;
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
offSet = 1;
Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');
inds = bwmorph(inds,'dilate',20);
Xf = Xf(:,inds);posd = posd(inds,:);veld = veld(inds,:);vel = vel(inds);pos = pos(inds,:);
Xf = bsxfun(@minus,Xf,mean(Xf,2));
Xf = zscore([real(Xf);imag(Xf)],0,2);
Xf = complex(Xf(1:end/2,:),Xf(end/2+1:end,:));
t = Xf;clear Xf;
%t = abs(t);%[real(t);imag(t)];
% %%FOR 1D TRACK
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
w = watershed(b==0);
w = w-1; 
posd(mod(w,2) ==1 ,1) = posd(mod(w,2) ==1 ,1) + max(posd(:));
pos(mod(w,2) ==1 ,1) = pos(mod(w,2) ==1 ,1) + max(pos(:));
runs1 = round(w/2);
inds = runs1 > 0 & runs1 <= max(runs);
t1 = zeros(size(t,1),max(runs),accumbins(1)*2);
for j = 1:size(t,1)
         t1(j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) 2*accumbins(1)] ,@mean);
end
spatial = randn(size(t1,1),2*accumbins(1));
for i = 1:size(t1,1)
    temp = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    [u,s,v] = svds(temp,1);
    spatial(i,:) = exp(1i*angle(mean(u)))*s*v';% 
%     if -min(v) > max(v) 
%         v = -v; 
%     end
end
if ~exist('posInds','var') || isempty(posInds)
    posInds = find(max(abs(spatial)') > 0);
end
t = t(posInds,:);t1 = t1(posInds,:,:);
t = reshape(zscore(t(:)),size(t));
t1 = reshape(zscore(t1(:)),size(t1));
figure;plot(abs(spatial)');
spatial = spatial(posInds,:);
[~,peakLoc] = max(abs(spatial)');
% [~,indLoc] = sort(peakLoc);
% %posInds = posInds(indLoc);
% spatial = spatial(indLoc,:);
% t = t(indLoc,:);t1 = t1(indLoc,:);
% peakLoc = peakLoc(indLoc);
meanAng = zeros(1,size(t,1));
for i = 1:size(t,1)
%    meanAng(i) = mean(r1(1:size(Xf,1)-1,posInds(i)));%
    meanAng(i) = mean(t(i,posd(:,1) == peakLoc(i)));% < peakLoc(i)+accumbins(1)/10 & posd(:,1) > peakLoc(i)-accumbins(1)/10));
    meanAng(i) = angle(meanAng(i));
%    meanAng(i) = circ_mean(angle(r1(1:size(Xf,1)-1,posInds(i))),abs(r1(1:size(Xf,1)-1,posInds(i))));
end
spatial = bsxfun(@times,spatial,exp(1i*(pi/2-meanAng')));
figure;imagesc(complexIm(spatial,0,1));
spsh = spatial;
for i = 1:size(spsh,1)
    spsh(i,:) = circshift(spsh(i,:),[0 -peakLoc(i)+accumbins(1)]);
end
figure;imagesc(complexIm(spsh,0,1));
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
tes = zeros(numel(posInds),max(runs),2*accumbins(1));
meanAng = zeros(1,size(t,1));
figure;
for i = 1:size(t,1)
    meanAng(i) = mean(t(i,posd(:,1) == peakLoc(i)));% < peakLoc(i)+accumbins(1)/10 & posd(:,1) > peakLoc(i)-accumbins(1)/10));
    meanAng(i) = angle(meanAng(i));
end
for i = 1:numel(posInds)
    te = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    subplot(xdim,ydim,i);imagesc(complexIm(imfilter(te.*exp(1i*(pi/2-meanAng(i))),fspecial('gaussian',5,1)),0,1,[],[],3));axis off;%s(temp(indLoc(i))),[0 max(te(:))]
    tes(i,:,:) = te;
    if exist('shank','var')
        title(shank(i));
    end
end

%sPlot([bsxfun(@times,t,sk'); vel']);
%figure;imagesc(complexIm(corr(ups(:,:)'),0,1));
superImpC(tes,[],1,prctile(abs(tes(:)),99));

% figure;
% for i = 1:max(runs1)
%     for j = 1:size(t,1)
%         indUse = runs1' == i & posd(:,1) < peakLoc(j)+accumbins(1)/10 & posd(:,1) > peakLoc(j)-accumbins(1)/10;
% %        plot(pos(indUse,1),posd(indUse,1));
% %        meanAng(j) = angle(mean(t(j,indUse)));
%         scatter(pos(indUse,1),angle(t(j,indUse)*exp(-1i*meanAng(j))),abs(t(j,indUse))*40,'filled');%angle(spatial(i,peakLoc(j)))
%  %       plot(pos(indUse,1),angle(t(j,indUse)*exp(-1i*meanAng(j))),'k');%angle(spatial(i,peakLoc(j)))
%         hold all;
%     end
%     hold off;pause(2);
% end