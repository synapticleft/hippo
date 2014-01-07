function [A,W,Z] = runTriggerICA(pos,Xf,thresh)
%% run complex-valued ICA; alternatively
%% convert demodulated complex data to real valued w/ 2x dimensionality, 
%% run fastICA
%warning off all;
dec = 1;
%bounds = [.1 .9];
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
for i = 1:size(pos,2)
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
if size(Xf,2) < size(pos,1)
    pos = pos(1:size(Xf,2),:);
end
nanInds = isnan(pos(:,1));
if size(pos,2) > 2
    nanInds = nanInds | isnan(pos(:,3));
end
pos = pos(~nanInds,:);Xf = Xf(:,~nanInds);%v = v(~nanInds,:);sp = sp(:,~nanInds);
if dec > 1
    for i = 1:4
        posd(:,i) = decimate(pos(:,i),dec);
    end
    pos = posd;clear posd;
end
% pos = bsxfun(@minus,pos,mean(pos));
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
%if dec > 1
%Xfd = decimate(Xf(1,:),dec);Xfd(2:size(Xf,1),:) = 0;
%for i = 2:size(Xf,1)
%    Xfd(i,:) = decimate(Xf(i,:),dec);
%end
%Xf = Xfd;clear Xfd;
%end

vel = angVel(pos);%vel = filtLow(vel(:,1),1250/32,1);
vel = [0; vel(:,1)];
vel = filtLow(vel,1250/32/dec,.5);
vel = vel/max(vel);
inds = vel > thresh;
%inds = abs(zscore(abs(v(:,1)))) < 1;
% sum(inds)
%Xf = Xf(:,inds);
%A = zeros(100,64,63);
%W = zeros(100,63,64);
%Z = zeros(100,64,63);
%for i = 1:100
%[A(i,:,:),W(i,:,:),Z(i,:,:)] = ACMNsym(Xf(:,inds),'mle_circ');%%%%nonCircComplexFastICAsym(Xf,'pow');%c%complex_ICA_EBM(Xf);%%zscore(Xf,0,2));%zscore(Xf,0,2)n
%end
%[A W Z] = cfpa2(Xf);%
[A W Z ] = ACMNsym(Xf(:,inds),'mle_circ');%alphascfastica(Xf(:,inds));%
%[W,sp] = grunica1(Xf(:,inds));
%W = W*sp; A = pinv(W);
% return
% Xf = [real(Xf);imag(Xf)];%[abs(Xf); angle(Xf)];
% rdim = size(Xf,1);
% [A,W,Z] = gfastica(zscore(Xf,0,2),'lastEig',rdim,'g','tanh','approach','symm','stabilization','on');%
%t = W*zscore(Xf,0,2);
