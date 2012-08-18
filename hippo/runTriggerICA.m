function [r t] = runTriggerICA(pos,v,Xf,accumbins,thresh,r)
%% convert demodulated complex data to real valued w/ 2x dimensionality, 
%% run fastICA, then bin and render activations.
warning off all;
dec = 1;
bounds = [.1 .9];
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
pos = bsxfun(@minus,pos,mean(pos));
%[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i))+.001);
    posd(:,i) = floor(pos(:,i)*accumbins(min(numel(accumbins),i)))+1;
end
%figure;imagesc(accumarray(posd,ones(size(posd,1),1),accumbins));
offSet = 1;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
     [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
for i = 2:size(Xf,1)
    Xfd(i,:) = decimate(Xf(i,:),dec);
end
Xf = Xfd;clear Xfd;
Xf = [real(Xf);imag(Xf)];%[abs(Xf); angle(Xf)];
%Xf = zscore(Xf,0,2);
%Xf = filtLow(Xf,1250/32/dec,2);
vel = filtLow(vel,1250/32/dec,1);
aV = decimate(abs(v(:,1)),dec);
aV = filtLow(aV,1250/32/dec,1);
aV = aV/prctile(aV,99.3);
%figure;plot(vel(:,1));hold all;plot(aV);return
%inds = aV > thresh & aV < 1;
vel = vel/max(vel);inds = vel > thresh;
%inds = inds & pos(:,1)< bounds(2) & pos(:,1) > bounds(1);
%sum(inds)
%t = Xf(:,inds);
%[r,~,t] = runica(Xf(:,inds),'pca',50);
rdim = size(Xf,1)-2;
if ~exist('r','var')
    %r = runica(zscore(Xf(:,inds),0,2),'pca',rdim);
    %[~,r] = gestimate(zscore(Xf(:,inds),0,2),rdim,500);%
    [~,r] = fastica(zscore(Xf(:,inds),0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');
end
r1 = pinv(r);
t = r*zscore(Xf,0,2);%(:,inds)
%figure;plot(abs(skewness(t')));drawnow;
sPlot([pos(inds,1:2)'; t(:,inds)],[],[],1);
xdim = ceil(sqrt(size(t,1)));ydim = ceil(size(t,1)/xdim);
% h1 = figure; h1s = figure;h2 = figure;
% for i = 1:size(t,1)
% figure(h1);subplot(xdim,ydim,i);imagesc(accumarray(posd(inds,1:2),t(i,inds),[],@mean)');axis off;
% figure(h1s);subplot(xdim,ydim,i);imagesc(accumarray(posd(inds,1:2),t(i,inds),[],@std)');axis off;
% figure(h2);subplot(xdim,ydim,i);subplot(xdim,ydim,i);imagesc(complexIm(reshape(...
%     complex(r1(1:32,i),r1(end-31:end,i)),[8 4]),0,1));axis off;
%     %complex(r1(1:64,i),r1(66:129,i)),[8 8]),0,1));axis off;
% end
% drawnow;
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
for k = 1:2
    runs1 = bwlabel(w>0 & mod(w,2) == k-1 & w <=2*max(runs));%b*((-1)^k)>0);
    inds = runs1 > 0;
    for j = 1:size(t,1)
        vInterp(k,j,:,:) = accumarray([runs1(inds); posd(inds,1)']',t(j,inds),[max(runs) accumbins(1)] ,@mean);
    end
end
t1 = [squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))];
h1 = figure;
for i = 1:size(t1,1)
    temp = reshape(t1(i,:),[max(runs) 2*accumbins(1)]);
    [~,s,v] = svds(temp,1);
    v = s*v'; 
    if -min(v) > max(v) 
        v = -v; 
    end
    spatial(i,:) = s*v';
    figure(h1);subplot(xdim,ydim,i);imagesc(temp);axis off;
%    figure(h2);subplot(xdim,ydim,i);imagesc(complexIm(reshape(complex(r1(1:32,i),r1(34:65,i)),[8 4]),0,1));axis off;
end
figure;plot(spatial');
temp = find(max(spatial') > 300);
[~,peakLoc] = max(spatial(temp,:)');
[~,indLoc] = sort(peakLoc);
figure;imagesc(spatial(temp(indLoc),:))