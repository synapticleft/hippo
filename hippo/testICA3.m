function [xc s] = testICA3(nChannels,SNRs,freq,nSources,trials,figsOn) 
%% simulation for complex valued data
div = 100;
pos = linspace(0,1,100);
v = exp(-1i*linspace(0,numel(pos),numel(pos)));
v = repmat(v,[1 trials]); v = v(:).';
sigma = [.5 1]*3;
[xs ys] = meshgrid(1:nChannels(2),1:nChannels(1));
Xf = zeros(prod(nChannels),trials*numel(pos));
if nSources > div
    scaleSources = floor(nSources/div);
    nSources = div;
else
    scaleSources = 1;
end
mix = zeros(nSources,prod(nChannels));
xcLen = round(numel(pos)/4);
acts = zeros(nSources,trials*numel(pos));
xc = zeros(4,xcLen*2+1);
sActs = 0;
for ii = 1:scaleSources
    tic;
    for i = 1:nSources
        randns = randn(2,trials);
        temp = exp(-((xs-rand*(nChannels(2)+1)).^2/sigma(1).^2 + (ys-rand*(nChannels(1)+1)).^2/sigma(2).^2));
        mix(i,:) = temp(:);
        r = rand;
        rs = max(rand/100,.0005);
        act = max(.1,1+randn/3)*giveCmor(pos'-r,freq,rs).';%max(0,filtLow(rand(numel(pos),1),numel(pos),10));
        temp1 = (max(.1,1+randns(1,:)*SNRs(1))'*act)';
        %for j = 1:trials
        %    acts(i,(j-1)*numel(pos)+(1:numel(pos))) = max(.1,1+randns(1,j)*SNRs(1))*circshift(act,[0 round(randns(2,j)*SNRs(2))]);
        %end
        acts(i,:) = temp1(:);
    end
    for j = 1:size(acts,1)
        xc(1,:) = xc(1,:) + xcov(acts(j,:),xcLen);
    end
     Xf = Xf + mix'*acts;
     sActs = sActs + sum(acts);
    toc
end
Xf = bsxfun(@times,Xf,v);
Xf = Xf + complex(randn(size(Xf)),randn(size(Xf)))*std(Xf(:))*SNRs(3);
[A,W,~] = ACMNsym(Xf,'mle_circ');%cfastica(Xf);%
Am = mean(A);
A = bsxfun(@times,A,exp(1i*-angle(Am)));
W = bsxfun(@times,W,exp(1i*-angle(Am)'));
Xf = bsxfun(@times,Xf,conj(v));
temp = W*Xf;
Xfs = squeeze(mean(reshape(Xf,[size(Xf,1) numel(pos) trials]),3));
%for i = 1:size(Xf,1)
%    xc(2,:) = xc(2,:) + xcov(Xfs(i,:),xcLen);
%end
temp = reshape(temp,[size(temp,1) numel(pos) trials]);
temp(:) = zscore(temp(:));
t = permute(temp,[1 3 2]);temp = squeeze(mean(t,2));
thresh = 2;
tInds = sum(abs(temp) > thresh,2) > numel(pos)/100;
[~,m] = max(abs(temp(tInds,:)),[],2);
%[~,m] = sort(m);
f = find(tInds);
%fn = find(~tInds);
%[~,mn] = max(abs(temp1(~tInds,:)),[],2);
%[~,mn] = sort(mn);
%A = A(:,[f(m); fn(mn)]);
%W = W([f(m);fn(mn)],:);
%A = A(:,1:numel(f));A = bsxfun(@minus,A,mean(A,2));
s = sum(abs(temp(tInds,:)));
acts = acts(1:25,:);
actsm = squeeze(mean(reshape(acts,[size(acts,1) numel(pos) trials]),3));
actss = std(reshape(acts,[size(acts,1) numel(pos) trials]),0,3);
actss = permute(reshape(actss,[size(acts,1),1,size(actss,2)]),[3 2 1]);
xs = repmat(1:numel(pos),[size(acts,1) 1]);
for i = 1:numel(f) 
    xc(3,:) = xc(3,:) + xcov(temp(f(i),:),xcLen);
end
ys = repmat(eye(numel(pos)),[1 trials]);
Xf = bsxfun(@minus,Xf,mean(Xf,2));
W = (Xf*Xf'+eye(size(Xf,1))*0)\Xf*ys';
ys = W'*Xf;
ys = squeeze(mean(reshape(ys,[size(ys,1) numel(pos) trials]),3));
for i = 1:size(ys,1)
    xc(4,:) = xc(4,:) + xcov(ys(i,:),xcLen);
end
xc = bsxfun(@rdivide,xc,max(abs(xc),[],2));
if figsOn
tot = 3;
figure;showGrid(t);
figure;subplot(1,tot,1);boundedline(xs',abs(actsm'),actss,'alpha');axis tight;
set(gca,'fontsize',16);title('Source activations');xlabel('Position');ylabel('activity');
%[Xz,Z] = zca2(Xf);
%Xzstd = std(reshape(Xz,[size(Xz,1) numel(pos) trials]),0,3);
%Xzs = Z*Xfs;elecs = [12 28 44];
%Xz = Xz(elecs,:);
%Xzs1 = interp1((1:numel(pos)),Xzs(elecs,:)',(1:.1:numel(pos)),'spline')';
%Xzstd1=  interp1(1:numel(pos),mean(Xzstd)',1:.1:numel(pos),'spline')';
%subplot(1,tot,2);
%scatter3(Xzs1(1,:),Xzs1(2,:),Xzs1(3,:),ones(1,size(Xzs1,2))*50,colormap(jet(size(Xzs1,2))),'filled');hold all;
col = colormap(jet(numel(pos)));
for i = 1:numel(f)
    cc(i,:) = col(m(i),:);
%    scatter3(Xzs(elecs(1),m),Xzs(elecs(2),m),Xzs(elecs(3),m),120,'k','filled');
end
%plot(zscore(real(Xzs(elecs(2),:))),'r','linewidth',2);hold all;plot(zscore(real(Xfs(elecs(2),:))),'b','linewidth',2);
%axis tight;
subplot(1,tot,2);
%set(gca,'nextPlot','add','ColorOrder',cc);
xs = repmat(1:numel(pos),[numel(f) 1]);
temps = std(t(tInds,:,:),0,2);
boundedline(xs',abs(temp(tInds,:)'),permute(temps,[3 2 1]),'cmap',cc);axis tight;
set(gca,'fontsize',16);title('ICA activations');xlabel('position');
subplot(1,tot,3);
plot(-xcLen:xcLen,abs(xc(1,:)),'b','linewidth',2);hold all;
%plot(abs(xc(2,:)),'b*');
plot(-xcLen:xcLen,abs(xc(3,:)),'r','linewidth',2);
plot(-xcLen:xcLen,abs(xc(4,:)),'r*','linewidth',2);axis tight;
legend({'Sources','ICA','Supervised'});
set(gca,'fontsize',16);xlabel('distance');ylabel('autocorrelation');title('Tuning')
%superImpC1(t,1,.05);
end

function mor = giveCmor(pos,Fc,Fb)
mor = exp(2*1i*pi*Fc*pos).*exp(-(pos.*pos)/Fb);