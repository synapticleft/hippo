function t = testICAFit(nChannels,SNRs,freq,nSources,trials) 
%% some older version of testICA
div = 100;
pos = linspace(0,1,100);% linspace(0,1,2000) linspace(0,1,4000)];
if ~freq
    v = ones(1,numel(pos));%real(exp(1i*linspace(0,numel(pos),numel(pos))));%
else
    v = exp(1i*linspace(0,numel(pos),numel(pos)));
end
v = repmat(v,[trials 1]); v = v(:).';
sigma = [.5 2]*3;
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
for ii = 1:scaleSources
    tic;
    for i = 1:nSources
        randns = randn(2,trials);
        temp = exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));
        mix(i,:) = temp(:);
        r = rand;
        rs(i) = max(rand/100,.0005);
        act = giveCmor(pos'-r,freq,rs(i)).';%exprnd(1)*
        for j = 1:trials
            acts(i,(j-1)*numel(pos)+(1:numel(pos))) = max(.1,1+randns(1,j)*SNRs(1))*circshift(act,[0 round(randns(2,j)*SNRs(2))]);
        end
    end
    %acts = acts + complex(randn(size(acts)),randn(size(acts)))*std(v)/SNRs(1);%sqrt(mean(v*v'))
    if ~freq
        acts = real(acts);
    end
    Xf = Xf + mix'*acts;
    toc
end
% xs = 0:.05:1;
% acts = reshape(acts,[size(acts,1) numel(pos) trials]);
% actsm = repmat(mean(acts,3),[1 1 trials]);
% acts = bsxfun(@rdivide,acts(:,:),sqrt(sum(acts(:,:).*conj(acts(:,:)),2)));
% actsm = bsxfun(@rdivide,actsm(:,:),sqrt(sum(actsm(:,:).*conj(actsm(:,:)),2)));
% figure;plot(xs,hist(abs(diag(acts*actsm')),xs));
% return
Xf = bsxfun(@times,Xf,v);
%Xf = bsxfun(@minus,Xf,mean(Xf,2));
if ~freq
    [a,b,c] = runica(Xf);
    figure;subplot(411);plot(abs(acts)');
    subplot(412);plot(abs(whiten(acts)'));
    subplot(413);plot(Xf');
    subplot(414);plot(c');
else
    Xf = Xf + complex(randn(size(Xf)),randn(size(Xf)))*std(Xf(:))*SNRs(3);
    t = linspace(0,1,1000);
    %%%tiled
    nBas = 5;
    cSt = 0;
    cEnd = 1;
    db = (cEnd - cSt) / (nBases-1)
    c = cSt:db:cEnd;
    bas = nan(length(phi),length(t));
    for k = 1:length(phi)
        bas(k,:) = ( cos( ...
            max(-pi, min(pi,pi*(t - c(k))/(db) )) ) ...
            + 1) / 2;
    end
    %%%
    [A,W,~] = cfastica(Xf);%ACMNSym(Xf,'mle_circ');%
    Am = mean(A);
    A = bsxfun(@times,A,exp(1i*-angle(Am)));
    A1 = bsxfun(@rdivide,A,sqrt(sum(A.*conj(A))));
    A1 = bsxfun(@minus,A1,mean(A1.').');
    figure;showGrid(A1,nChannels);
    Xf = bsxfun(@times,Xf,conj(v));
    temp = W*Xf;
    temp = reshape(temp,[size(temp,1) numel(pos) trials]);t = permute(temp,[1 3 2]);temp = squeeze(mean(t,2));
    thresh = 1;
    tInds = sum(abs(temp) > thresh,2) > numel(pos)/100;
    [~,m] = max(abs(temp(tInds,:)),[],2);
    [~,m] = sort(m);
    f = find(tInds);
    cc = abs(corr(temp([f(m); find(~tInds)],:)'));
    for i = 1:numel(f)
        ccD(i) = mean(diag(cc(1:numel(f),1:numel(f)),i));
    end
    figure;subplot(211);imagesc(cc);subplot(212);
    figure;subplot(311);plot(abs(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3))');
    subplot(312);plot(abs(mean(reshape(whiten(Xf),[size(Xf,1) numel(pos) trials]),3)'));axis tight;
    subplot(313);plot(abs(temp(tInds,:)'));axis tight;%hold all;plot(abs(xcorr
    superImpC1(t,1,.05);
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
    figure;imagesc(hist3([xs(:) allMs(:)],[100 100]));
    figure;hist(xs(:)-allMs(:),100);
    %subplot(414);hold all;
    %for i = 1:numel(tInds)
    %   absInds = abs(temp(tInds(i),:)) > thresh / 2;
    %   scatter(find(absInds),angle(temp(tInds(i),absInds)),abs(temp(tInds(i),absInds))*20,'filled');
    %end
end

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*