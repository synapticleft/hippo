function medErr = testICA3r(shanks,SNRs,smoothness,nSources,trials) %[xc fs]
%% simulation for paper using real-valued, non-oscillatory activity

figsOn = 1;
div = 50;
pos = linspace(0,1,500);
posHet = pos;%[linspace(0,.66,100) linspace(.67,1,100)];
nChannels = shanks;
sigma = 2;
Xf = zeros(prod(nChannels),trials*numel(pos));
if nSources > div
    scaleSources = floor(nSources/div);
    nSources = div;
else
    scaleSources = 1;
end
mix = zeros(nSources,prod(nChannels));
acts = zeros(nSources,trials*numel(pos));
XfIm = zeros(size(Xf,1),numel(pos),3);
rgbCols = ones(size(acts,1),1,3);
gWinPos = gausswin(smoothness);%50 default;gWinElec = gausswin(10);
xcLen = round(numel(pos)/4);
xc = zeros(3,xcLen*2+1);

for ii = 1:scaleSources
    tic;
    for i = 1:nSources
        randns = randn(2,trials);
        %temp = randn(nChannels);%exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));%randn(nChannels);%
        temp = exp(-(linspace(1,prod(nChannels),prod(nChannels))-(rand*1.1-.05)*prod(nChannels)).^2/sigma(1).^2);
        %temp = max(0,filter(gWinElec,1,randn(round(prod(nChannels)*1.2),1)));
        mix(i,:) = temp(:);%(round(.15*prod(nChannels))+(1:prod(nChannels)));%
        rs = max(rand/100,.0005);
        %act = giveCmor(pos'-rand,freq,rs).';%
        %act = max(-0,filtLow(randn(numel(pos)*1.2,1),numel(pos)*1.2,10))';%max(.1,1+randn/3)*%max(0,);%zeros(1,numel(pos));act(ceil(rand*numel(pos))) = 1;%randn(1,numel(pos));%exprnd(1)*
        act = max(0,filter(gWinPos,1,randn(numel(pos)*1.2,1)))';
        act = act(:,floor(.15*numel(pos))+(1:numel(pos)));
        act = max(0,interp1(pos,act,posHet));
        temp1 = (max(.1,1+randns(1,:)*SNRs(1))'*act)';
        acts(i,:) = temp1(:);
    end
    actsM = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3));
    for i = 1:size(acts,1)
        xc(1,:) = xc(1,:) + xcov(actsM(i,:),xcLen);
    end
     Xf = Xf + mix'*acts;%
     actsMean = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3));
     [~,m] = max(actsMean');
     rgbCols(:,:,1) = m/numel(pos);
     cols = [1 1 1];%squeeze(hsv2rgb(rgbCols));%
     for i = 1:3
        XfIm(:,:,i) = XfIm(:,:,i) + bsxfun(@times,mix,cols(:,i))'*actsMean;
    end
    toc
end
[~,s] = sort(m);
actsIm = makeIm(actsMean,m/numel(pos),s);
mixIm = makeIm(mix,m/numel(pos),s);
XfIm = min(1,XfIm/prctile(XfIm(:),99.9));
if figsOn
figure;subplot(131);image(permute(mixIm,[2 1 3]));
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel 'neuron #';ylabel 'electrode #';
subplot(132);image(actsIm);
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel 'position';ylabel 'neuron #';
subplot(133);image(XfIm);
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel 'position';ylabel 'electrode #';
end
Xf = bsxfun(@minus,Xf,mean(Xf,2));
r = randn(size(Xf));%filtLow(,numel(pos),3);
Xf = Xf + r/std(r(:))*std(Xf(:))*SNRs(2);
ys = repmat(eye(numel(pos)),[1 trials]);
%Xf = bsxfun(@minus,Xf,mean(Xf,2));
tr = rand(1,numel(pos)*trials) < .5;
te = ~tr;
[u,s] = eig(Xf*Xf');
poses = repmat(linspace(1,numel(pos),numel(pos))',[trials 1]);
if 1
Xfa = pinv(sqrt(s))*u'*Xf;
for i = 1:size(Xf,1)
    Xf = Xfa(end-i+1:end,:);
W = (Xf(:,tr)*Xf(:,tr)'+eye(size(Xf,1)))\Xf(:,tr)*ys(:,tr)';
[~,m] = max( W'*Xf(:,te));
medErr(i) = median(abs(m'-poses(te)));
end
else
    W = (Xf(:,tr)*Xf(:,tr)'+eye(size(Xf,1)))\Xf(:,tr)*ys(:,tr)';
    [~,m] = max( W'*Xf(:,te));
end
return
ys = W'*Xf;
ys = squeeze(mean(reshape(ys,[size(ys,1) numel(pos) trials]),3));
for i = 1:size(ys,1)
    xc(3,:) = xc(3,:) + xcov(ys(i,:),xcLen);
end
if figsOn
ctrs{1} = 1:numel(pos);ctrs{2} = ctrs{1};
temp = hist3([m' poses(te)],ctrs);
temp = temp(:);
[xs ys] = meshgrid(1:numel(pos),1:numel(pos));
figure;scatter(xs(temp >0),ys(temp>0),temp(temp>0),'k','s','filled');axis square;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);ylabel('predicted');xlabel('actual');
end
[A W] = fastica(Xf,'approach','symm');%ACMNsym(Xf,'mle_circ');%
%[A W] = ACMNsym(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))),'mle_circ');
%[A W] = fastica(bsxfun(@times,Xf,real(exp(1i*linspace(pi,pi*numel(pos)*trials,numel(pos)*trials)))),'approach','symm');
sk = sign(skewness((W*(Xf))')');
c = bsxfun(@times,W*(Xf),sk);
A = bsxfun(@times,A,sk');
%timeCourses = c(:,1:numel(pos)*2);
c = reshape(c,[size(c,1) numel(pos) trials]);%b*Xf;
c = permute(c,[1 3 2]);
temp = c;temp(:,:) = zscore(temp(:,:),0,2);%temp(:) = zscore(temp(:));
c = squeeze(mean(temp,2));
[mx] = max(c,[],2);
f = find(mx > 2);% & mp ~= 1 & mp ~= numel(pos));%max(c(:))/2);
fs = numel(f);
%nf = find(max(c,[],2) <=2);
%[~,snf] = sort(m1(nf));
%figure;set(gca,'nextPlot','add','ColorOrder',cc);
%timeCourses = timeCourses(f(sf),:);
%plot(bsxfun(@plus,timeCourses,linspace(0,-size(timeCourses,1)*3,size(timeCourses,1))')','linewidth',2);
if ~figsOn
figure;
    col = colormap(hsv(numel(pos)));
    for i = 1:size(c,1)
        [~,m] = max(c(i,:));
        cc(i,:) = col(m,:);
    end
    %figure;showGrid(temp(f,:,:));
set(gca,'nextPlot','add','ColorOrder',cc);
%plot(c','linewidth',2);axis tight;
xs = repmat(posHet,[numel(f) 1]);
temps = std(temp(f,:,:),0,2);
subplot(211);boundedline(xs',(c(f,:)'),permute(temps,[3 2 1]),'cmap',cc(f,:));axis tight;
subplot(212);boundedline(xs',bsxfun(@rdivide,c(f,:),max(c(f,:),[],2))',permute(temps,[3 2 1])/10,'cmap',cc(f,:));axis tight;
set(gca,'fontsize',16);xlabel 'position'; ylabel 'activation';drawnow;
end
ca = c;
for i = 1:size(ca,1)
    ca(i,:) = interp1(posHet,ca(i,:),pos);
end
if figsOn
    [~,m1] = max(c');[~,sf] = sort(m1(f));
figure;subplot(121);image(permute(makeIm((A)',m1/numel(pos),f(sf)),[2 1 3]));
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel 'source #';ylabel 'electrode #';
subplot(122);image(makeIm((ca),m1/numel(pos),f(sf)));%axis image;
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel 'position';ylabel 'source #';%; nf(snf)
end
for i = 1:numel(f) 
    xc(2,:) = xc(2,:) + xcov(c(f(i),:),xcLen);
end
xc = bsxfun(@rdivide,xc,max(abs(xc),[],2));
if figsOn
figure;
plot(-xcLen:xcLen,xc(1,:),'b','linewidth',2);hold all;
plot(-xcLen:xcLen,xc(2,:),'r','linewidth',2);
plot(-xcLen:xcLen,xc(3,:),'r*','linewidth',2);axis tight;
legend({'Sources','Unsupervised','Supervised'});
set(gca,'fontsize',16,'xtick',[],'ytick',[]);xlabel('distance');ylabel('autocorrelation');title('Tuning');
end

function actsIm = makeIm(acts,col,ord)
actsIm= zeros(size(acts,1),size(acts,2),3);
acts = max(0,acts - prctile(acts(:),1));
actsIm(:,:,3) = 1-min(1,acts/prctile(acts(:),99.9));
actsIm(:,:,2) = 0;%1;
for i = 1:size(acts,1)
    actsIm(i,:,1) = col(i);
end
actsIm = actsIm(ord,:,:);
actsIm(:,:,3) = flipud(squeeze(actsIm(:,:,3)));
actsIm = hsv2rgb(actsIm);

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);%((pi*Fb)^(-0.5))*