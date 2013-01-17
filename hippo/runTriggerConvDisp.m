function [ac numSamps] = runTriggerConvDisp(X,pos,phi,thresh,accumbins,dewhiteningMatrix)

ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
peakToPeak = ceil(1250/dec/8);
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
vela = vel;
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
%vel = vel/max(vel);
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
if 1
    if ~exist('dewhiteningMatrix','var')
        numSamples = 50000;
        indsSub = rand(numel(inds),1) < numSamples/sum(inds);
        [Ex, Dx] = eig(cov(X(:,indsSub' & inds)'));
        d = flipud(diag(Dx));
        cumVar = sum(d);
        maxLastEig = sum(cumsum(d)/cumVar < .9999999)
        Dx = Dx(end-maxLastEig+1:end,end-maxLastEig+1:end);
        Ex = Ex(:,end-maxLastEig+1:end);
        factors = diag(Dx);
        noise_factors = ones(size(Dx,1),1);
        rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .999999)
        noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind)));
        Dx = diag(factors./noise_factors);
        whiteningMatrix = sqrt(inv(Dx)) * Ex';
        dewhiteningMatrix = Ex * sqrt (Dx);
    else
        whiteningMatrix = pinv(dewhiteningMatrix);
    end
        X = whiteningMatrix * X;
%          phi1 = phi;
%     [~,J,R] = size(phi);
%     phi = zeros(size(whiteningMatrix,1),J,R);
%     for j = 1:J
%         phi(:,j,:) =  whiteningMatrix * squeeze(phi1(:,j,:));
%     end
else
    X = bsxfun(@rdivide,X,std(X,0,2));
end

opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 0.01, 'cb', @cb_a);
[N,J,R] = size(phi);
ac =0;numSamps = 0;
lambda = [.5 3/50];%2.5;
%la = linspace(-3,2,10);%[.1 .5 1 2 5];
%laa = linspace(-5,2,10);%[.005 .01 .02 .05 1];
%la = linspace(la(5),la(7),10);
%laa = linspace(laa(6),laa(8),10);
%la = linspace(-.5,.2,8);
%laa = linspace(-1.1,.1,8);
la = linspace(.1,1,8);
laa = linspace(.005,.05,8);
%la = .5;laa = 3/50;
la = .8; laa = .25;
figure;
for k = 1:numel(la)
    for l = 1:numel(laa)
for j = 40%1:max(reg)
Xsamp = X(:,reg == j);
    %% compute the map estimate
    tic
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    a0 = zeros(J, P);
    %% no bounds
    lb  = zeros(1,J*P); % lower bound
    ub  = zeros(1,J*P); % upper bound
    nb  = ones(1,J*P); % bound type (none)
    [a1,fx,exitflag,userdata] = lbfgs(@objfun_a_conv, a0(:), lb, ub, nb, opts_lbfgs_a, Xsamp, phi, [la(k) laa(l)]);
    a1 = reshape(a1, J, P);
    %[~,id] = meshgrid(1:S,1:J);id = id';
    %aTemp = a1(:,1:S);aTemp = aTemp';
    %rpos = repmat(posd(reg == j),[J 1]);
    %plot(id(:));hold all;plot(rpos);plot(aTemp(:));return
    %[size(id) size(id(:)) size(rpos) size(rpos(:))]
    %ac = ac + accumarray([id(:) rpos],aTemp(:),[J max(posd)],@sum);
    %numSamps = numSamps + accumarray(posd(reg == j),ones(1,sum(reg == j)),[max(posd) 1],@sum);
    %sPlot(a1,[],0);
subplot(numel(la),numel(laa),(k-1)*numel(laa)+l);imagesc(a1);axis off tight;
    %subplot(211);imagesc(a1);
    %subplot(212);imagesc(bsxfun(@rdivide,ac,numSamps'));
    drawnow;
    time_inf = toc;

%     %% reconstruct
%     EI = zeros(N,S);
%     for r = 1:R
%         EIr = phi(:,:,r)*a1;
%         srt = R+1-r;
%         fin = R+S-r;
%         EI = EI + EIr(:,srt:fin);
%     end
end
    end
end
