function [A,W,Z] = runTriggerICA(pos,v,Xf,thresh)
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
for i = 1:size(pos,2)
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
%plot(pos(:,1),'r');return;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:size(pos,2)
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
if size(pos,2) > 2
    nanInds = nanInds | isnan(pos(:,3));
end
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);%sp = sp(:,~nanInds);
if dec > 1
    for i = 1:4
        posd(:,i) = decimate(pos(:,i),dec);
    end
    pos = posd;clear posd;
end
pos = bsxfun(@minus,pos,mean(pos));
% [a,~,~] = svd(pos(:,1:2),'econ');pos = a;
% for i = 1:2    
%     pos(:,i) = pos(:,i) - min(pos(:,i));
%     pos(:,i) = pos(:,i)/(max(pos(:,i)));
%     pos(:,i) = min(pos(:,i),.9999);
% end
% offSet = 1;
%Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).')];
%Xf = filtLow(Xf,1250/32,2);
%Xf = bsxfun(@times,Xf,exp(1i*angle(v(:,1)))');
if dec > 1
Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
for i = 2:size(Xf,1)
    Xfd(i,:) = decimate(Xf(i,:),dec);
end
Xf = Xfd;clear Xfd;
end

vel = angVel(pos);%vel = filtLow(vel(:,1),1250/32,1);
vel = [0; vel(:,1)];
vel = filtLow(vel,1250/32/dec,1);
vel = vel/max(vel);
inds = vel > thresh;
%inds = abs(zscore(abs(v(:,1)))) < 1;
sum(inds)
Xf = Xf(:,inds);
[A,W,Z] = ACMNsym(Xf,'mle_circ');%cfpa2(Xf);%%%nonCircComplexFastICAsym(Xf,'pow');%cfastica(Xf);c%complex_ICA_EBM(Xf);%%zscore(Xf,0,2));%zscore(Xf,0,2)n
return
 Xf = [real(Xf);imag(Xf)];%[abs(Xf); angle(Xf)];
 rdim = size(Xf,1);
 [A,W,Z] = gfastica(zscore(Xf,0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');%
%t = W*zscore(Xf,0,2);
