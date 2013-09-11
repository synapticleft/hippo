function [ac numSamps SAll] = runTriggerConvDisp1(X,pos,phi,thresh,accumbins,zp)
%% make plots of responses by convolutional sparse coding after learning
%% wrote for spatioTEMPORALLY whitened data


ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
%%Processing of position information
bounds = [.1 .9];
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
pos = pos(:,1);
for i = 1:size(pos,2)   
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
pos(nanInds) = 0;
vel = filtLow(vel,1250/32,.5);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
vel = vel/max(vel);
vel = resample(vel,ratio,1);
pos = resample(pos,ratio,1);
pos = pos(1:size(X,2),:); 
posd = floor(pos*accumbins(1))+1;posd = min(accumbins,max(1,posd));
vel = vel(1:size(X,2));
inds = vel > thresh(1);
%inds(pos < bounds(1) | pos > bounds(2)) = 0;
reg = bwlabel(inds);
h = hist(reg,0:max(reg));
a = accumarray(reg'+1,pos,[],@mean);
f = find(h(2:end) < 1250/dec*thresh(2) | a(2:end)' < bounds(1) | a(2:end)' > bounds(2));
inds(ismember(reg,f)) = 0;
reg = bwlabel(inds);
dPos = [0; diff(pos)];
a = accumarray(reg'+1,dPos,[],@mean);
f = find(a(2:end) > 0);
posd(ismember(reg,f)) = posd(ismember(reg,f)) + accumbins;
%% whiten X
X = bsxfun(@minus,X,mean(X,2));
f = find(inds);
[N,J,R] = size(phi);
Y = zeros(N,numel(f));
for i = 1:R
        Y = Y + zp(R/2*size(X,1)+(1:size(X,1)),(i-1)*N+(1:N))*X(:,f-R/2+i);
    %Y = Y + whiteningMatrix(:,(i-1)*N+(1:N))*X(:,f-R/2+i);%circshift(X,[0 i+R/2]);
end
X = Y;
reg = reg(f);posd = posd(f);
clear Y;
opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 0.01, 'cb', @cb_a);
ac =0;numSamps = 0;
%la = .2;laa = 3/50;
la = 1.2; laa = nan;
params.Fs = 1250/dec;params.tapers = [3 5];
figure;
SAll = 0;
for j = 1:max(reg)
Xsamp = X(:,reg == j);
    %% compute the map estimate
    tic
     S = size(Xsamp,2);
%     P = S+R-1;	% number of selection locations for source generation
%     a0 = zeros(J, P);
%     %% no bounds
%     lb  = zeros(1,J*P); % lower bound
%     ub  = zeros(1,J*P); % upper bound
%     nb  = ones(1,J*P); % bound type (none)
%     [a1] = lbfgs(@objfun_a_conv, a0(:), lb, ub, nb, opts_lbfgs_a, Xsamp, phi, [la laa]);% [la(k) laa(l)]);
%     a1 = reshape(a1, J, P);
     [~,id] = meshgrid(1:S,1:J);id = id';
%     aTemp = a1(:,1:S);aTemp = aTemp';
    aTemp = Xsamp';
    [S,f1] = mtspectrumc(aTemp,params);
    SAll = SAll + interp1(f1,S,0:.5:100);
    rpos = repmat(posd(reg == j),[J 1]);
    ac = ac + accumarray([id(:) rpos],abs(aTemp(:)),[J max(posd)],@sum);
    numSamps = numSamps + accumarray(posd(reg == j),ones(1,sum(reg == j)),[max(posd) 1],@sum);
    subplot(211);imagesc(bsxfun(@rdivide,ac,numSamps'));drawnow;
    subplot(212);imagesc(log(SAll));drawnow;
    drawnow;
    time_inf = toc;
end