function t = testSparse(shanks,SNRs,nSources,trials) 
%% attempt at overcomplete sparse coding of simulated data -- didnt work so far
div = 100;
pos = linspace(0,1,64);% linspace(0,1,2000) linspace(0,1,4000)];
nChannels = shanks;
sigma = [.5 1]*5;
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
randPos = sort(rand(1,nSources*scaleSources),'ascend');
for ii = 1:scaleSources
    tic;
    for i = 1:nSources
        randns = randn(2,trials);
        temp = exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));
        mix(i,:) = temp(:);
        r = randPos((ii-1)*nSources+i);%i/nSources;%rand;
        rs(i) = max(rand/100,.0005);
        act = randn*abs(giveCmor(pos'-r,0,rs(i))).';%exprnd(1)*
        for j = 1:trials
            acts(i,(j-1)*numel(pos)+(1:numel(pos))) = max(.1,1+randns(1,j)*SNRs(1))*circshift(act,[0 round(randns(2,j)*SNRs(2))]);
        end
    end
    Xf = Xf + mix'*acts;%
    toc
end
    Xf = Xf + randn(size(Xf))*std(Xf(:))*SNRs(3);
    %Xf = (bsxfun(@times,Xf,sin(linspace(0,size(Xf,2),size(Xf,2)))));
    %Xf = bsxfun(@minus,Xf,mean(Xf,2));
    %[A,W] = gfastica(Xf,'approach','symm');
    Xfa = Xf;
    [Xf,wh] = whiten(bsxfun(@times,Xf,real(exp(1i*linspace(pi,numel(pos)*trials*pi,numel(pos)*trials)))));
    de = pinv(wh);
    %Xf = bsxfun(@rdivide,Xf,std(Xf,0,2));
    %wh = eye(size(Xf,1));de = wh;
    %[Xf,wh,de] = whiten(Xf);
    unittest;
    W = pinv(phi); A = phi;
    figure;showGrid(de*A,nChannels);
    temp = abs(W*wh*Xfa);
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
    %figure;subplot(311);plot(abs(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3))');axis tight;
    %subplot(312);plot(abs(A'));axis tight;%mean(reshape(whiten(Xf),[size(Xf,1) numel(pos) trials]),3)'));axis tight;
    %subplot(313);plot(abs(temp(tInds,:)'));axis tight;%hold all;plot(abs(xcorr
    superImpC1(t,1,.05);

    
function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*