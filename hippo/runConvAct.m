function [magAll timeAll numAll numSamps] = runConvAct(X,pos,phi1)

load('phi1.mat','c','d');
spatial = zeros(size(c,1),size(c,3));
c2 = c;
for i = 1:size(c,1)
c2(i,:,:) = squeeze(c(i,:,:))./max(d,1);
[u,s,v1] = svds(squeeze(c2(i,:,:)),1);
spatial(i,:) = mean(u)*s*v1';
end
posInds = find(max(spatial') > .3*max(spatial(:)));
spatial = spatial(posInds,:);
[~,peakLoc] = max(abs(spatial)');
[~,indLoc] = sort(peakLoc);
posInds = posInds(indLoc);
spatial = spatial(indLoc,:);
%plot(spatial');
%superImp(c2,posInds,1);
%% position processing
thresh = [.05 1];accumbins = 50;
ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
%%Processing of position information
bounds = [.2 .9];
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
vel = [0; vel];
pos = bsxfun(@minus,pos,nanmean(pos));
[~,~,c] = svd(pos(~nanInds,1:2),'econ');pos = (c\pos(:,1:2)')';%pos = a;pos(nanInds) = nan;
pos = pos(:,1);
pos = pos - min(pos) + eps;
pos = pos/(max(pos)+eps);
pos(nanInds) = 0;
vel = filtLow(vel,1250/32,.5);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
%vel = vel(1:size(X,2));
inds = vel > thresh(1);
%% which repetition of rat running
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
%runs = bwlabel(b > 0);
runs = watershed(b==0);
f = find(runs == 0);
runs(f) = runs(f-1);
% for i = 1:max(runs)
%     if mean(diff(pos(w == i))) < 0
%         pos(w == i) = 2 - pos(w == i);
%     end
% end
%% which contiguous chunk of data
chunk = bwlabel(inds);
h = hist(chunk,0:max(chunk));
a = accumarray([ones(size(chunk)); chunk+1]',pos,[],@mean);
f = find(h(2:end) < 1250/dec*thresh(2) | a(2:end) < bounds(1) | a(2:end) > bounds(2));
inds(ismember(chunk,f)) = 0;
%chunk = bwlabel(inds);
chunk1 = bwlabel(round(resample(double(inds),ratio,1)));chunk1 = chunk1(1:size(X,2));
dPos = [0; diff(pos)];
a = accumarray(runs'+1,dPos,[],@mean);
f = find(a(2:end) > 0);
pos(ismember(runs,f)) = 2-pos(ismember(runs,f));pos = pos/2;
%posd = floor(pos*accumbins*2)+1;posd = min(2*accumbins,max(1,posd));
pos = resample(pos,ratio,1);
pos = pos(1:size(X,2)); 
posd1 = floor(pos*accumbins*2)+1;posd1 = min(accumbins*2,max(1,posd1));
runs = ceil(runs/2);
runs1 = round(resample([runs runs(end)*ones(1,100)],ratio,1));runs1 = runs1(1:size(X,2));
%% 1 basis function
folder = dir('../*J=001*');
load(['../' folder(1).name '/params.mat'],'whiteningMatrix','lambda');
load(['../' folder(1).name '/phi.mat'],'phi');
wh1 = whiteningMatrix;lambda1 = lambda;phi1 = phi;%phi1 = repmat(phi,[1 1 3]);
%% many bases
load('params.mat','lambda','whiteningMatrix','dewhiteningMatrix');
load('phi.mat','phi');
opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20,'factr', 0.01, 'cb', @cb_a);
[N,J,R] = size(phi);J1 = size(phi1,2);Jp = numel(posInds)
% phiA = zeros(size(dewhiteningMatrix,1),J,R);
% for j = 1:J
%     %phiA(:,j,:) = filtLow(squeeze(phi(:,j,:)),1250/8,12,2);%
%     phiA(:,j,:) = dewhiteningMatrix*squeeze(phi(:,j,:));
% end
% phi = phiA;
la = linspace(.1,1,8);
laa = linspace(.005,.05,8);
%la = la(1);
%laa = laa(3);
la = .6;laa = .25;%linspace(.05,.45,9);%.25;
magAll = zeros(2,Jp,2*accumbins);timeAll = magAll;numAll = magAll; numSamps = 0;
maxTime = 1250/dec/6/2;mag = .5;
figure;
for k = 1:numel(la)
    for l = 1:numel(laa)
for j = 1:max(chunk1)
    %% for single
    Xsamp = wh1*X(:,chunk1 == j);%
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    %% no bounds
%      a0 = zeros(J1, P);
%      a1 = lbfgs(@objfun_a_conv, a0(:), zeros(1,J1*P), zeros(1,J1*P), ones(1,J1*P), opts_lbfgs_a, Xsamp, phi1, lambda1);%[.2 .25]);%
%      a1 = a1';
    %% for multi
    Xsamp = whiteningMatrix*X(:,chunk1 == j);%
    xComp = angle(morFilter(Xsamp(end,:),8,1250/dec,1,.2));
    %Xsamp = filtLow(Xsamp,1250/8,12,2);
    lb  = zeros(1,J*P); % lower bound
    ub  = zeros(1,J*P); % upper bound
    nb  = 0*ones(1,J*P); % bound type (none)
    a0 = zeros(size(phi,2), P);
    aMult = lbfgs(@objfun_a_conv, a0(:), lb, ub, nb, opts_lbfgs_a, Xsamp, phi,lambda);% [la(k) laa(l)]);%
    aMult = reshape(aMult, J, P);
    %% display
    %subplot(numel(la),numel(laa),(k-1)*numel(laa)+l);
    %sPlot([aMult;a1*5],[],0);drawnow;
    %subplot(311);plot(X(:,chunk1 == j)');axis tight;
    %subplot(412);plot(reCon(phi,aMult,S,dewhiteningMatrix)');axis tight;
    aMult = aMult(posInds,:);
    %subplot(312);plot(reCon(phi(:,posInds,:),aMult,S,dewhiteningMatrix)');axis tight;
    %aMult = filtfilt(gausswin(3),1,aMult')';
    subplot(313);plot(aMult');axis tight;%imagesc(aMult,[0 prctile(aMult(:),99.9)]);axis off;colorbar;hold all;drawnow;%
%    plot(a1,'r','linewidth',2);
%    peakTimes = bwdist(a1 > mag);
%    dTimes = [diff(peakTimes) -1];dTimes(dTimes == 0) = -1;
%    peakTimes = peakTimes.*dTimes;peakTimes = max(-maxTime,min(maxTime,peakTimes));
%    peakTimes = repmat(peakTimes(1:S),[Jp 1])';
    peakTimes = repmat(exp(1i*xComp),[Jp 1])';
    [~,id] = meshgrid(1:S,1:Jp);id = id';
    aMult = aMult(:,1:S)';f = find(aMult(:) > mag);
    posTemp = repmat(posd1(chunk1 == j),[Jp 1]);
    magAll(1,:,:) = squeeze(magAll(1,:,:)) + accumarray([id(f) posTemp(f)],aMult(f),[Jp max(posd1)],@sum);
    timeAll(1,:,:) = squeeze(timeAll(1,:,:)) + accumarray([id(f) posTemp(f)],peakTimes(f),[Jp max(posd1)],@sum);
    numAll(1,:,:) = squeeze(numAll(1,:,:)) + accumarray([id(f) posTemp(f)],ones(numel(f),1),[Jp max(posd1)],@sum);
    f = find(aMult(:) < -mag);
    magAll(2,:,:) = squeeze(magAll(2,:,:)) + accumarray([id(f) posTemp(f)],aMult(f),[Jp max(posd1)],@sum);
    timeAll(2,:,:) = squeeze(timeAll(2,:,:)) + accumarray([id(f) posTemp(f)],peakTimes(f),[Jp max(posd1)],@sum);
    numAll(2,:,:) = squeeze(numAll(2,:,:)) + accumarray([id(f) posTemp(f)],ones(numel(f),1),[Jp max(posd1)],@sum);
    numSamps = numSamps + accumarray([ones(numel(posTemp),1) posTemp(:)],ones(numel(posTemp),1),[1 max(posd1)],@sum);
    subplot(311);imagesc([squeeze(magAll(1,:,:)) squeeze(magAll(2,:,:))]);
    subplot(312);imagesc(complexIm([squeeze(numAll(1,:,:).*exp(1i*angle(timeAll(1,:,:)))) squeeze(numAll(2,:,:).*exp(1i*angle(timeAll(2,:,:))))],0,1));drawnow; 
    %subplot(312);imagesc(complexIm([squeeze(timeAll(1,:,:)./max(1,numAll(1,:,:))) squeeze(timeAll(2,:,:)./max(1,numAll(2,:,:)))],0,1));drawnow;%input('');
%    ac1 = ac1 + accumarray([id(f) repmat(runs1(chunk1 == j)',[J 1]) ...
%        repmat(posd1(chunk1 == j),[J 1])],aMult(f),[J max(runs1) max(posd)],@sum);
 %   numSamps1 = numSamps1 + accumarray([runs1(chunk1 == j)' posd1(chunk1 == j)],ones(1,sum(chunk1 == j)),[max(runs1) max(posd)],@sum);
end
    end
end

function EI = reCon(phi,a1,S,de)
[N,J,R] = size(phi);
    EI = zeros(N,S);
    for r = 1:R
        EIr = phi(:,:,r)*a1;
        srt = R+1-r;
        fin = R+S-r;
        EI = EI + EIr(:,srt:fin);
    end
    if exist('de','var')
        EI = de*EI;
    end