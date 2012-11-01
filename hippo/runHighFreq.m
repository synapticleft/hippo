function [A,tes,td] = runHighFreq(X,pos,A,W)

accumbins = [50 1];
ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
%%Processing of position information
thresh = .05;bounds = [.1 .9];win = [-1 1]*ceil(1250/dec/8/2);
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
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,nanmean(pos));
[~,~,c] = svd(pos(~nanInds,1:2),'econ');pos = (c\pos(:,1:2)')';%pos = a;pos(nanInds) = nan;
for i = 1:size(pos,2)   
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
%pos(nanInds,:) = 0;
%%THE filtLow function requires some functions from the signal processing
%%toolbox, but is not particularly necessary.
vel = filtLow(vel,1250/32,1);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
vel = vel/max(vel);
vel = resample(vel,ratio,1);
pos = resample(pos,ratio,1);
pos = pos(1:size(X,2),:); vel = vel(1:size(X,2));
inds = vel > thresh;
pos = pos(inds,:);
for i= 1:2
    posd(:,i) = max(1,floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1);
end
if ~exist('A','var')
    [A,W] = gfastica(zscore(X(:,inds),0,2),'lastEig',size(X,1),'g','tanh','approach','symm','stabilization','on');
end
%%
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
w = watershed(b==0);
w = w-1; %w(w== max(w)) = 0;
t = W*zscore(X,0,2);%
td = zeros(size(X,1),ceil(size(X,2)/ratio));
for i = 1:size(X,1)
    td(i,:) = decimate(t(i,:).^2,ratio);
end
t = t(:,inds);
vInterp = zeros(2,size(t,1),max(runs),accumbins(1));
for k = 1:2
%    runs1 = b*(-1^k)>0;runs1 = bwlabel(runs1);
    runs1 = bwlabel(w>0 & mod(w,2) == k-1 & w <=2*max(runs));
    inds = runs1 > 0;
    for j = 1:size(t,1)
        vInterp(k,j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) accumbins(1)] ,@std);
    end
end
t1 = [squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))];
posInds = 1:size(t1,1);
xdim = ceil(sqrt(numel(posInds)));ydim = ceil(numel(posInds)/xdim);
%h1 = figure;h2 = figure;
tes = zeros(size(t1,1),max(runs),2*accumbins(1));
for i =1:size(t1,1)
    te = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    te(te == 0) = prctile(te(:),20);
    tes(i,:,:) = te;
%    figure(h1);subplot(xdim,ydim,i);imagesc(imfilter(te,fspecial('gaussian',5,1),'replicate'));axis off;
%    figure(h2);subplot(xdim,ydim,i);imagesc(reshape(A(:,i)*sign(skewness(A(:,i))),[8 numel(A(:,i))/8]));axis off;
end