function allResps = testICABatch(trials) 
%% early version when i was figuring out parameter space of ICA on simulations
freq = 10;
div = 100;
pos = linspace(0,1,100);
if ~freq
    v = ones(1,numel(pos));
else
    v = exp(1i*linspace(0,numel(pos),numel(pos)));
end
v = repmat(v,[trials 1]); v = v(:).';
nChannels = [8 8];
sigma = [.5 2]*3;
[xs ys] = meshgrid(1:nChannels(2),1:nChannels(1));
Xf = zeros(prod(nChannels),trials*numel(pos));
nSources = 100;%[10 100 1000 10000];
chanSNRs = [0 .01 .1 1];%[inf 100 10 1];
sourceSNRs = [0 1 2 4];
sourceSNRs1 = [0 1 2 4];
allResps = zeros(numel(nSources),numel(sourceSNRs),numel(sourceSNRs1),numel(chanSNRs),prod(nChannels),numel(pos));
for k = 1:numel(nSources)
    for kk = 1:numel(sourceSNRs)
        for kk1 = 1:numel(sourceSNRs1)
        Xf(:) = 0;
    if nSources(k) > div
        nSource = div;
        scaleSources = floor(nSources(k)/div);
    else
        nSource = nSources(k);
        scaleSources = 1;
    end
    mix = zeros(nSource,prod(nChannels));
    acts = zeros(nSource,trials*numel(pos));
    for ii = 1:scaleSources
        tic;
        for i = 1:nSource
            randns = randn(2,trials);
            temp = exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));
            mix(i,:) = temp(:);
            r = rand;
            rs(i) = max(rand/100,.0005);
            act = giveCmor(pos'-r,freq,rs(i)).';%exprnd(1)*
            for j = 1:trials
                acts(i,(j-1)*numel(pos)+(1:numel(pos))) = max(.1,1+randns(1,j)*sourceSNRs(kk))*circshift(act,[0 round(randns(2,j)*sourceSNRs1(kk1))]);
            end
        end
        %acts = acts + complex(randn(size(acts)),randn(size(acts)))*std(v)/SNRs(1);%sqrt(mean(v*v'))
        if ~freq
            acts = real(acts);
        end
        Xf = Xf + mix'*acts;
        toc
    end
    Xf = bsxfun(@times,Xf,v);
    [V,D] = eig(Xf*Xf'/size(Xf,2));
    for kkk = 1:numel(chanSNRs)
        if 1
        [~,W] = cfastica(Xf + V*sqrt(D)*complex(randn(size(Xf)),randn(size(Xf)))*chanSNRs(kkk));%ACMNSym(Xf,'mle_circ');%*std(Xf(:))
        temp = W*bsxfun(@times,Xf,conj(v));
        temp = reshape(temp,[size(temp,1) numel(pos) trials]);t = permute(temp,[1 3 2]);
        allResps(k,kk,kk1,kkk,:,:) = squeeze(mean(t,2));
        else
            acts = reshape(acts,[size(acts,1) numel(pos) trials]);
            actsm = repmat(mean(acts,3),[1 1 trials]);
            acts = bsxfun(@rdivide,acts(:,:),sqrt(sum(acts(:,:).*conj(acts(:,:)),2)));
            actsm = bsxfun(@rdivide,actsm(:,:),sqrt(sum(actsm(:,:).*conj(actsm(:,:)),2)));
            allResps1(kk,kk1,kkk) = mean(abs(diag(acts*actsm')));
        end
    end
        end
    end
end
allResps = squeeze(allResps);

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*