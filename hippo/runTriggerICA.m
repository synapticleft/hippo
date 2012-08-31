function [W t] = runTriggerICA(pos,v,Xf,accumbins,thresh)
%% convert demodulated complex data to real valued w/ 2x dimensionality, 
%% run fastICA, then bin and render activations.
warning off all;
dec = 1;
bounds = [.1 .9];
pos(pos == -1) = nan;
%figure;plot(pos(:,1));hold all;
%plot([0; diff(pos(:,1))]);
%plot(flipud([0; diff(flipud(pos(:,1)))]));
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
if dec > 1
    for i = 1:4
        posd(:,i) = decimate(pos(:,i),dec);
    end
    pos = posd;clear posd;
end
vel = angVel(pos);%vel = filtLow(vel(:,1),1250/32,1);
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
if dec > 1
Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
for i = 2:size(Xf,1)
    Xfd(i,:) = decimate(Xf(i,:),dec);
end
Xf = Xfd;clear Xfd;
end
Xf = [real(Xf);imag(Xf)];%[abs(Xf); angle(Xf)];
%Xf = zscore(Xf,0,2);
%Xf = filtLow(Xf,1250/32/dec,2);
vel = filtLow(vel,1250/32/dec,1);
%aV = decimate(abs(v(:,1)),dec);
%aV = filtLow(aV,1250/32/dec,1);
%aV = aV/prctile(aV,99.3);
%figure;plot(vel(:,1));hold all;plot(aV);return
%inds = aV > thresh & aV < 1;
vel = vel/max(vel);inds = vel > thresh;
%inds = inds & pos(:,1)< bounds(2) & pos(:,1) > bounds(1);
%sum(inds)
%t = Xf(:,inds);
%[r,~,t] = runica(Xf(:,inds),'pca',50);
rdim = size(Xf,1)/2;
[~,W] = fastica(zscore(Xf(:,inds),0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');%
t = W*zscore(Xf,0,2);
