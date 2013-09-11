function t = testICA2(shanks,SNRs,freq,nSources,trials) 
%% older clunky version of simulation
div = 100;
pos = linspace(0,1,200);% linspace(0,1,2000) linspace(0,1,4000)];
if ~freq
    v = ones(1,numel(pos));%real(exp(1i*linspace(0,numel(pos),numel(pos))));%
else
    v = exp(1i*linspace(0,numel(pos),numel(pos)));
end
v = repmat(v,[1 trials]); v = v(:).';
nChannels = shanks;
sigma = [.5 1]*2;
[xs ys] = meshgrid(1:nChannels(2),1:nChannels(1));
Xf = zeros(prod(nChannels),trials*numel(pos));
if nSources > div
    scaleSources = floor(nSources/div);
    nSources = div;
else
    scaleSources = 1;
end
mix = zeros(nSources,prod(nChannels));
acts = zeros(nSources,trials*numel(pos));
randPos = rand(1,nSources*scaleSources);%sort(rand(1,nSources*scaleSources),'ascend');
xc = zeros(4,numel(pos)+1);
for ii = 1:scaleSources
    tic;
    for i = 1:nSources
        randns = randn(2,trials);
        %complex(randn(nChannels),randn(nChannels));%% 
        %temp = randn(nChannels);%exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));%randn(nChannels);%
        temp = exp(-(linspace(1,prod(nChannels),prod(nChannels))-((ii-1)*nSources+i)/(nSources*scaleSources)*prod(nChannels)).^2/(max(.1,1+randn*sigma(1))).^2);
        mix(i,:) = temp(:);
        r = rand;%i/nSources;%randPos((ii-1)*nSources+i);%
        rs(i) = max(rand/100,.0005);
        act = max(-0,filtLow(randn(numel(pos)*1.2,1),numel(pos)*1.2,10))';%max(.1,1+randn/3)*giveCmor(pos'-r,freq,rs(i)).';%max(0,);%zeros(1,numel(pos));act(ceil(rand*numel(pos))) = 1;%randn(1,numel(pos));%exprnd(1)*
        act = act(:,.1*numel(pos)+(1:numel(pos)));
        temp1 = (max(.1,1+randns(1,:)*SNRs(1))'*act)';
        acts(i,:) = temp1(:);
%         for j = 1:trials
%             acts(i,(j-1)*numel(pos)+(1:numel(pos))) = max(.1,1+randns(1,j)*SNRs(1))*circshift(act,[0 round(randns(2,j)*SNRs(2))]);
%         end
    end
%     for j = 1:size(acts,1)
%         xc(1,:) = xc(1,:) + xcov(acts(j,:),numel(pos)/2);
% %         for jj = 1:size(acts,1)
% %             if j ~= jj
% %                xc(2,:) = xc(2,:) + xcov(acts(j,:),acts(jj,:),numel(pos)/2);
% %             else
% %                 xc(3,:) = xc(3,:) + xcov(acts(j,:),acts(jj,:),numel(pos)/2);
% %             end
% %         end
%     end
     Xf = Xf + mix'*acts;%
    toc
end
% for j = 1:size(Xf,1)
%     xc(4,:) = xc(4,:) + xcov(Xf(j,:),numel(pos)/2);
% end
% xc(4,:) = xc(4,:)/max(xc(4,:))*max(xc(1,:));
% figure;plot(xc');return
if ~freq
    Xf = real(Xf);
    acts = real(acts);
end
% xs = 0:.05:1;
% acts = reshape(acts,[size(acts,1) numel(pos) trials]);
% actsm = repmat(mean(acts,3),[1 1 trials]);
% acts = bsxfun(@rdivide,acts(:,:),sqrt(sum(acts(:,:).*conj(acts(:,:)),2)));
% actsm = bsxfun(@rdivide,actsm(:,:),sqrt(sum(actsm(:,:).*conj(actsm(:,:)),2)));
% figure;plot(xs,hist(abs(diag(acts*actsm')),xs));
% return
Xf = bsxfun(@times,Xf,v);
Xf = bsxfun(@minus,Xf,mean(Xf,2));
if ~freq
    fica = 1;
    r = randn(size(Xf));%filtLow(,numel(pos),3);
    Xf = Xf + r/std(r(:))*std(Xf(:))*SNRs(3);
    if fica
        %[~,wh] = zca2(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))));
        %wh = eye(size(Xf,1));
    [A W] = fastica(Xf,'approach','symm');%wh*bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))));sp = wh;
    %[A W] = fastica(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))),'approach','symm');
    %[W,sp] = runica(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))),'sphering','on');%);%
    %[W,sp] = grunica1(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))));
    %W = W*sp;
    else
    %Xf = bsxfun(@minus,Xf,mean(Xf,2));
    params.n=size(Xf,1);   % dimensionality of input
%[Xf1, mean_patch, wh] = preprocess(bsxfun(@times,Xf,real(exp(1i*linspace(0,numel(pos)*trials,numel(pos)*trials)))));
%[Xf1,wh,de] = whiten(bsxfun(@times,Xf,real(exp(1i*linspace(pi,numel(pos)*trials*pi,numel(pos)*trials)))),0,0);
%
Xf1 = bsxfun(@times,Xf,real(exp(1i*linspace(0,numel(pos)*trials,numel(pos)*trials))));
Xf1 = bsxfun(@rdivide,Xf1,std(Xf1,0,2));
figure;imagesc(Xf1*Xf1');drawnow;
wh = eye(size(Xf,1));
%m = sqrt(sum(Xf.^2) + (1e-8));10
%Xf = bsxfunwrap(@rdivide,Xf,m);
%% Run the optimization
params.lambda = .1;
params.numFeatures = 100;
params.epsilon = 1e-5;
%configure minFunc
options.Method = 'lbfgs';
options.MaxFunEvals = Inf;
options.MaxIter = 300;
%options.display = 'off';
%options.outputFcn = 'showBases';
% initialize with random weights
randTheta = randn(params.numFeatures,params.n)*0.01;  % 1/sqrt(params.n);
randTheta = randTheta ./ repmat(sqrt(sum(randTheta.^2,2)), 1, size(randTheta,2)); 
randTheta = randTheta(:);
% optimize
[opttheta, cost, exitflag] = minFunc( @(theta) softICACost(theta, Xf1/10, params), randTheta, options);   % Use x or xw 
W = reshape(opttheta, params.numFeatures, params.n);
%    [c,b,a] = gfastica(Xf,'approach','symm');
%    Xf = squeeze(mean(reshape(Xf,[size(Xf,1) numel(pos) trials]),3));
    end
    sk = skewness((W*(Xf))')';
    c = bsxfun(@times,W*(Xf),sk);
    c = reshape(c,[size(c,1) numel(pos) trials]);%b*Xf;
    c = permute(c,[1 3 2]);
    figure;showGrid(c);
    figure;subplot(2,1,1);plot(squeeze(std(c,0,2))');subplot(2,1,2);plot(squeeze(mean(c,2))');
    %c1 = c;
    c = squeeze(mean(c,2));
    c(:) = zscore(c(:));
    Xfs = squeeze(mean(reshape(Xf,[size(Xf,1) numel(pos) trials]),3));
    ca = bsxfun(@times,W*(Xfs),sk);%pinv(sp)*
    figure;subplot(1,5,1);showGrid(mix(1:4,:),nChannels);
    subplot(1,5,2);plot(squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3))','linewidth',2);axis tight;
    Xfs1 = interp1((1:numel(pos)),Xfs(1:3,:)',(1:.1:numel(pos)))';
    subplot(1,5,3);scatter3(Xfs1(1,:),Xfs1(2,:),Xfs1(3,:),ones(size(Xfs1,2),1)*20,colormap(jet(size(Xfs1,2))),'filled');hold all;
    col = colormap(jet(numel(pos)));
    for i = 1:size(c,1)
        [v,m] = max(c(i,:));
        cc(i,:) = col(m,:);
        if v > 2
            scatter3(Xfs(1,m),Xfs(2,m),Xfs(3,m),120,cc(i,:),'filled');
            scatter3(Xfs(1,m),Xfs(2,m),Xfs(3,m),120,'k','linewidth',2);
        end
    end
    subplot(1,5,4);
    set(gca,'nextPlot','add','ColorOrder',cc);
    plot(c','linewidth',2);axis tight;
    f = find(max(c,[],2) > 2);%max(c(:))/2);
    c = c(f,:);ca = ca(f,:);
    %figure;plot(squeeze(std(c1(f,:,:),0,2))');
    for i= 1:numel(f)
        [~,m] = max(c(i,:));
        c(i,:) = circshift(c(i,:),[0 numel(pos)/2-m]);
        ca(i,:) = circshift(ca(i,:),[0 numel(pos)/2-m]);
    end
    [~,z] = zca2(Xf(:,:));
    Xfsw = z*Xfs;
    for i = 1:size(Xf,1)
        c1(i,:) = xcov(Xfs(i,:),numel(pos)/2,'coeff');
        c1a(i,:) = xcorr(Xfs(i,:),numel(pos)/2,'unbiased');
        c1b(i,:) = xcorr(Xfsw(i,:),numel(pos)/2,'coeff');
    end
    acts = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3));
    for i = 1:size(acts,1)
        c2(i,:) = xcov(acts(i,:),numel(pos)/2,'coeff');
    end
    cm = mean(c);c1a = mean(c1a(:,2:end));c1a = c1a/max(c1a);
    c1b = mean(c1b(:,2:end));
    subplot(1,5,5);plot(cm,'r','linewidth',2);hold all;
    plot(c1b*max(cm),'r--','linewidth',2);
    %plot(mean(ca)/max(mean(ca))*max(cm),'r--','linewidth',2);
    plot(mean(c1(:,2:end))*max(cm),'k','linewidth',2);
    plot(c1a*max(cm),'k--','linewidth',2);
    plot(mean(c2(:,2:end))*max(cm),'g','linewidth',2);axis tight;
    %ys = (Xfs'*Xfs)\Xfs'*Xf(:,:);
    ys = repmat(eye(numel(pos)),[1 trials]);
    Xf = bsxfun(@minus,Xf,mean(Xf,2));
    W = (Xf*Xf'+eye(size(Xf,1)))\Xf*ys';
    ys = W'*Xf;
    ys = squeeze(mean(reshape(ys,[size(ys,1) numel(pos) trials]),3));
    for i = 1:size(ys,1)
        ys(i,:) = circshift(ys(i,:),[0 size(ys,2)/2-i]);
    end
    plot(mean(ys)/max(mean(ys))*max(cm),'b--','linewidth',2);
    %mx = max(c') < max(-c');
    %c(mx,:) = -c(mx,:);
    %subplot(411);plot(abs(acts)');
    %subplot(412);plot(abs(whiten(acts)'));
    %subplot(413);plot(Xf');
    %subplot(414);plot(c');
else
    %[x,whiteningMatrix,dewhiteningMatrix] = whiten(squeeze(mean(reshape(bsxfun(@times,Xf,conj(v)),[size(Xf,1) numel(pos) trials]),3)));
    %D = eig(Xf*Xf'/size(Xf,2));
    Xf = Xf + complex(randn(size(Xf)),randn(size(Xf)))*std(Xf(:))*SNRs(3);
    Xfs = squeeze(mean(reshape(Xf,[size(Xf,1) numel(pos) trials]),3));
    Xfs = repmat(Xfs,[1 trials]);
    [U,D] = svd(Xf,'econ');
    D1 = abs(U'*(Xf*Xf')*U);
    D2 = abs(U'*(Xfs*Xfs')*U); 
    %D1 = flipud(eig(Xf*Xf'/size(Xf,2)));
    %D2 = flipud(eig(Xfm*Xfm'/size(Xfm,2)));
    Xfn = Xf - Xfs;
    D3 = abs(U'*(Xfn*Xfn')*U);
    %D3 = flipud(eig(Xfm*Xfm'/size(Xfm,2)));
    %figure;%plot(log(max(eps,D)));
    %hold all;plot(log(max(eps,diag(D.^2))));plot(log(max(eps,diag(D1))));hold all;plot(log(max(eps,diag(D2))));plot(log(max(eps,diag(D3))));plot(log(1+diag(D2)./diag(D3)));
    %[~,whiteningMatrix,dewhiteningMatrix] = whiten(Xf1);
    %Xf2 = whiteningMatrix*Xf;
    [A,W,~] = cfastica(Xf);%ACMNsym(Xf,'mle_circ');%
    %[W,sp] = grunica1(Xf,'on');W = W*sp;A = pinv(W);
    %A = dewhiteningMatrix * W';  
  	%W = W * whiteningMatrix;
    Am = mean(A);
    A = bsxfun(@times,A,exp(1i*-angle(Am)));
    A1 = bsxfun(@rdivide,A,sqrt(sum(A.*conj(A))));A = A1;
    %A1 = bsxfun(@minus,A1,mean(A1.').');
    %figure;showGrid(A,nChannels);
    %[u,~] = eig(Xf*Xf');
    %v = (u(:,end)\Xf);
    Xf = bsxfun(@times,Xf,conj(v));
    temp = W*Xf;
    temp = reshape(temp,[size(temp,1) numel(pos) trials]);
    temp(:) = zscore(temp(:));
    t = permute(temp,[1 3 2]);temp = squeeze(mean(t,2));
    temp1 = squeeze(mean(abs(t),2));
    thresh = 1.5;
    tInds = sum(abs(temp) > thresh,2) > numel(pos)/100;
    [~,m] = max(abs(temp(tInds,:)),[],2);
    [~,m] = sort(m);
    f = find(tInds);
    fn = find(~tInds);
    [~,mn] = max(abs(temp1(~tInds,:)),[],2);
    [~,mn] = sort(mn);
    A = A(:,[f(m); fn(mn)]);
    W = W([f(m);fn(mn)],:);
    %figure;imagesc(abs(W*(Xf*Xf')*W'));
    D1 = max(eps,abs(diag(W*(Xfs*Xfs')*W')));
    D2 = max(eps,abs(diag(W*(Xfn*Xfn')*W')));
    %plot(log(D1));plot(log(D2));plot(log(1+D1./D2));
    figure;showGrid((t([f(m);fn(mn)],:,:)));
    A = A(:,1:numel(f));A = bsxfun(@minus,A,mean(A,2));
    figure;showGrid(A',shanks);
    %figure;plot(kurtosis(real(t([f(m);fn(mn)],:))'));hold all;
    %plot(kurtosis(imag(t([f(m);fn(mn)],:))'));plot(kurtosis(abs(t([f(m);fn(mn)],:))'));
    figure;subplot(221);imagesc(abs(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3)));
    subplot(222);imagesc(abs(A));
    subplot(223);imagesc(abs(temp([f(m); fn(mn)],:)));
    subplot(224);imagesc(temp1([f(m); fn(mn)],:));
    cc = abs(corr(temp([f(m); fn(mn)],:)'));
    for i = 1:numel(f)
        ccD(i) = mean(diag(cc(1:numel(f),1:numel(f)),i-1));
    end
    %figure;subplot(211);imagesc(cc);subplot(212);plot(ccD)
    figure;subplot(311);plot(abs(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3))');axis tight;
    subplot(312);plot(abs(A'));axis tight;%mean(reshape(whiten(Xf),[size(Xf,1) numel(pos) trials]),3)'));axis tight;
    subplot(313);plot(abs(temp(tInds,:)'));axis tight;%hold all;plot(abs(xcorr
    superImpC1(t,1,.05);
    acts = reshape(acts,[size(acts,1) numel(pos) trials]);acts = permute(acts,[1 3 2]);
    tm = abs(squeeze(mean(t(f(m),:,:),2)));
    ts = squeeze(std(t(f(m),:,:),0,2));
    am = abs(squeeze(mean(acts,2)));
    as = squeeze(std(acts,0,2));
    figure;subplot(221);plot(tm');subplot(222);plot(ts');subplot(223);plot(am');subplot(224);plot(as');
    figure;scatter(tm(:),ts(:));hold all;scatter(am(:),as(:));
    Xf = reshape(Xf,[size(Xf,1) numel(pos) trials]);
    Xfm = squeeze(mean(Xf,3));
    Xfm = bsxfun(@rdivide,Xfm,sqrt(sum(Xfm.*conj(Xfm))));
    %figure;
    allMs = zeros(trials,numel(pos));
    xs = meshgrid(1:numel(pos),1:trials);
    for i = 1:trials
        m = abs(Xfm'*squeeze(Xf(:,:,i)));
        %subplot(211);imagesc(m);
        [~,m1] = max(m);
        %subplot(212);scatter(1:numel(pos),m1,'filled');hold all;drawnow;
        allMs(i,:) = m1;
    end
    %figure;imagesc(hist3([xs(:) allMs(:)],[numel(pos) numel(pos)]));
    %figure;hist(xs(:)-allMs(:),100);
    %subplot(414);hold all;
    %for i = 1:numel(tInds)
    %   absInds = abs(temp(tInds(i),:)) > thresh / 2;
    %   scatter(find(absInds),angle(temp(tInds(i),absInds)),abs(temp(tInds(i),absInds))*20,'filled');
    %end
end

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*