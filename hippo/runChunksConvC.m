function wFilt1 = runChunksConvC(X,v,pos,act,A,wFilt,wh)
% act and A need to be phase-centered.

ratio = round(size(X,2)/size(pos,1));
dec = 32/ratio;
peakToPeak = ceil(1250/dec/8);
%%Processing of position information
thresh = [.1 1];bounds = [.1 .9];win = [-1 1]*ceil(1250/dec/8);
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
%%THE filtLow function requires some functions from the signal processing
%%toolbox, but is not particularly necessary.
vel = filtLow(vel,1250/32,.5);
vel = vel/max(vel);
vel = resample(vel,ratio,1);pos = resample(pos,ratio,1);

pos = pos(1:size(X,2),:); vel = vel(1:size(X,2));
inds = vel > thresh(1);
reg = bwlabel(inds);
h = hist(reg,0:max(reg));
a = accumarray(reg+1,pos,[],@mean);
f = find(h(2:end) < 1250/dec*thresh(2) | a(2:end)' < bounds(1) | a(2:end)' > bounds(2));
inds(ismember(reg,f)) = 0;
reg = bwlabel(inds);
X = wh*X;
if ~exist('wFilt','var') || isempty(wFilt)
    act = resample(act.',ratio,1).';act = [zeros(size(act,1),size(X,2)-size(act,2)) act];
    wFilt = zeros(size(X,1),size(act,1),peakToPeak);
%wFiltI = wFilt;h = 0;h1 = h;
xCo = zeros(size(act,1),peakToPeak,peakToPeak);%xCoI = xCo;
for i = 1:max(reg)
    runInds = find(reg == i);
%        oldTheta = angle(v(runInds));%*exp(1i*-1*pi/6));
%        oldTheta = unwrap(oldTheta);
%        newTheta = linspace(oldTheta(1),oldTheta(end),round((oldTheta(end)-oldTheta(1))/2/pi*peakToPeak));
%        d = [0 diff(mod(newTheta,2*pi))];
%        d = find(d < -pi);
%        newTheta = newTheta(d(1):d(end)-1);
%        newY = interp1(oldTheta,X(:,runInds)',newTheta)';
%        newX = interp1(oldTheta,act(:,runInds).',newTheta).';
    for j = 1:size(act,1)
        t = [0 diff(angle(act(j,runInds)))] < -pi;
        f = [find(t) numel(t)];t = double(t);
        for k = 1:numel(f)-1
            t(f(k)) = max(0,mean(abs(act(j,runInds(1)-1+(f(k):f(k+1)))))-2);
        end
        t = toeplitz(zeros(1,peakToPeak),t);
        wFilt(:,j,:) = squeeze(wFilt(:,j,:)) + X(:,runInds)*t';
        xCo(j,:,:) = squeeze(xCo(j,:,:)) + t*t';
%         t = [0 diff(angle(newX(j,:)))] < -pi;
%         f = [find(t) numel(t)];t = double(t);
%         for k = 1:numel(f)-1
%             t(f(k)) = max(0,mean(abs(newX(j,f(k):f(k+1))))-2);
%         end
%         t = toeplitz(zeros(1,peakToPeak),t);%newX(j,:));
%         wFiltI(j,:,:) = squeeze(wFiltI(j,:,:)) + newY*t';
%         xCoI(j,:,:) = squeeze(xCoI(j,:,:)) + t*t';
    end
end
for i = 1:size(wFilt,2)
   wFilt(:,i,:) = (squeeze(xCo(i,:,:))\squeeze(wFilt(:,i,:))')';
%   wFilt(:,i,:) = sqrt(sum(sum(wFilt(:,i,:).^2)));
%   wFiltI(i,:,:) = squeeze(wFiltI(i,:,:))/squeeze(xCoI(i,:,:));
end
else
    dFilt = wFilt;wFilt = zeros(size(wh,1),size(wFilt,2),size(wFilt,3));
    for j = 1:size(wFilt,2)
        wFilt(:,j,:) =  wh * squeeze(dFilt(:,j,:));
        wFilt(:,j,:) = wFilt(:,j,:)/sqrt(sum(sum(wFilt(:,j,:).^2)));
    end
end
f1 = figure;
f2 = figure;
de = pinv(wh);
showPhi(dFilt,A,f1);
%for j = 1:50
peakToPeak = size(dFilt,3);
wFilt1 = zeros(size(X,1),size(wFilt,2),peakToPeak);
xCo = zeros(size(wFilt,2),peakToPeak,peakToPeak);
figure;
x = [.5 1 2 3];y = [1 3 5 7 9]/50;
%for a = 1:4
%    for b = 1:5
%        subplot(4,5,(a-1)*5+b);
for i = 1:max(reg)
    runInds = find(reg == i);
    t1 = inferC(X(:,runInds),wFilt,[.5 3/50]);%,[x(a) y(b)]
%    wFilt = updatePhi(X(:,runInds),wFilt,t1);
    t1 = t1(:,1:numel(runInds));
%    figure(f2);imagesc(t1);colorbar;drawnow;
   for j = 1:size(wFilt1,2)
       t = toeplitz([t1(j,1) zeros(1,peakToPeak-1)],t1(j,:));
       wFilt1(:,j,:) = squeeze(wFilt1(:,j,:)) + X(:,runInds)*t';
       xCo(j,:,:) = squeeze(xCo(j,:,:)) + t*t';
   end
%    showPhi(wFilt,A,f1);
end
%    end
%end
%end
reg = mean(xCo(:))*10;
for i = 1:size(wFilt,2)
   wFilt1(:,i,:) = ((squeeze(xCo(i,:,:)) + reg*eye(size(xCo,3)))\squeeze(wFilt1(:,i,:))')';
end
    dFilt = wFilt1;wFilt1 = zeros(size(de,1),size(wFilt1,2),size(wFilt1,3));
    for j = 1:size(wFilt1,2)
        wFilt1(:,j,:) =  de * squeeze(dFilt(:,j,:));
        wFilt1(:,j,:) = wFilt1(:,j,:)/sqrt(sum(sum(wFilt1(:,j,:).^2)));
    end
showPhi(wFilt1,A,f2);

function showPhi(wFilt,A,f)
figure(f);
xdim = ceil(sqrt(size(A,2)));ydim = ceil(size(A,2)/xdim);
for i = 1:size(wFilt,2)%size(Wp,2)-1
    subplot(xdim,ydim,i);
    set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(A(:,i),0,1)));
    plot(squeeze(wFilt(:,i,:))');axis tight;
%     figure(f2);subplot(xdim,ydim,i);
%     set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(A(:,i),0,1)));
%     plot(squeeze(wFiltI(i,:,:))');axis tight;
%     [S,f,Sc] = mtspectrumc(temp,params);
%     figure(f3);subplot(xdim,ydim,i);set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(A(:,i),0,1)));
%     plot(f,S);axis tight;drawnow;%;set(gca,'xtick',[],'ytick',[]);axis tight;%title(numel(highResp));drawnow;%cc(min(size(W,2),i))
end

function phi = updatePhi(Xsamp,phi,a1)
[N,J,R] = size(phi);eta = .01;
            [obj0,g] = objfun_phi(phi(:), Xsamp, a1);
            dphi = reshape(g, N, J, R);
            phi = phi - eta * dphi;

function a1 = inferC(Xsamp,phi,lambda)
if ~exist('lambda','var')
    lambda = [2 .5];%[2 .5] awesome;
end
opts_lbfgs_a = lbfgs_options('iprint', -1, 'maxits', 20, 'factr', 0.01, 'cb', @cb_a);%);
    [~,J,R] = size(phi);
    %% compute the map estimate
    tic
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    a0 = zeros(J, P);
    %% no bounds
    lb  = zeros(1,J*P); % lower bound
    ub  = zeros(1,J*P); % upper bound
    nb  = zeros(1,J*P); % bound type (none)
    
    a1 = lbfgs(@objfun_a_conv, a0(:), lb, ub, nb, opts_lbfgs_a, Xsamp, phi, lambda);
    a1 = reshape(a1, J, P);
    imagesc(a1);drawnow;