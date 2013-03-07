%% panel 3 - Demodulation demo
nSteps = 1000;
x = linspace(0,20*pi,nSteps);
y = exp(1i*x);
A = [zeros(200,1); ones(300,1); zeros(500,1)]/2+1;
win = 20;
A = filtfilt(gausswin(win),sum(gausswin(win)),A);
dP = [zeros(600,1); ones(300,1); zeros(100,1)]*-pi/2;
dP = filtfilt(gausswin(win),sum(gausswin(win)),dP);
z = A'.*exp(1i*(x+dP'));zd = -conj((z.*conj(y))*1i);
angCol = colormap('hsv');
c = ceil((angle(z)+pi)/(2*pi)*64);
figure;plot3(1:nSteps,real(y),imag(y),'k');
hold all;scatter3(1:nSteps,real(z),imag(z),[],bsxfun(@times,angCol(c,:),(A/max(A))),'filled');
sh = 4;
plot3(1:nSteps,real(y),imag(y)-sh,'k');
c = max(1,ceil((angle(zd)+pi)/(2*pi)*64));
hold all;scatter3(1:nSteps,real(zd),imag(zd)-sh,[],bsxfun(@times,angCol(c,:),(A/max(A))),'filled');
set(gca,'xtick',[],'ytick',[],'ztick',[],'linewidth',2);
%% panel 4/5
inds = 52965:53140;angCol = colormap('hsv');
%temp = filtLow(angVel(pos)',1250/32,2);
[posa,s,u] = svds(pos(:,1:2),1);
posa = s*posa;
temp = filtLow(diff(posa(inds)),1250/32,2);
indsa = inds(1)*4:inds(end)*4;
X1 = morFilter(X(:,indsa(1)-1000:indsa(end)+1000),8,1250/8);X1 = X1(:,1001:end-1000);
[u,s,v1] = svds(X1,1);
um = mean(abs(u)).*exp(1i*circ_mean(angle(u)));
v2 = um*s*v1';
%% panel 4
sub = 100:285;mA = max(max(abs(X1(:,sub))));
v1c = mean(abs(v2(sub)))*exp(1i*angle(v2(sub)));
figure;plot3((sub-min(sub))/1250*8,real(v1c),imag(v1c),'k');hold all;
sh = s/50;
plot3((sub-min(sub))/1250*8,real(v1c),imag(v1c)-sh,'k');hold all;
%X1d = bsxfun(@times,X1(:,sub),exp(1i*(-angle(v1c))));
X1d = exp(1i*angle(conj(u*s*v1(sub)'))).*X1(:,sub);
X1d = abs(X1d).*exp(1i*(angle(X1d) + pi/2));
for i = 1:size(X1,1)
    c = ceil((angle(X1(i,sub))+pi)/(2*pi)*64);
    scatter3((sub-min(sub))/1250*8,real(X1(i,sub)),imag(X1(i,sub)),[],bsxfun(@times,angCol(c,:),abs(X1(i,sub)')/mA),'filled');
    c = min(64,max(1,ceil((angle(X1d(i,:))+pi)/(2*pi)*64)));
    scatter3((sub-min(sub))/1250*8,real(X1d(i,:)),imag(X1d(i,:))-sh,[],bsxfun(@times,angCol(c,:),abs(X1(i,sub)')/mA),'filled');
end
%% panel 5
figure;subplot(411);plot(temp*1250/32,'k','linewidth',2);axis tight;
colorbar;set(gca,'xtick',[],'fontsize',16);title('Velocity');ylabel('cm/s');
subplot(412);imagesc(X(:,indsa));set(gca,'ytick',[1 64],'xtick',[],'fontsize',16);colormap jet;colorbar;freezeColors;title('Raw LFP');
subplot(413);imagesc(complexIm(X1,0,1));colorbar;set(gca,'ytick',[1 64],'xtick',[],'fontsize',16);ylabel('Channel #');title('Filtered Theta');
subplot(414);imagesc(linspace(0,size(X1,2)/1250*8,size(X1,2)),1:64,...bsxfun(@times,X1,exp(1i*angle(v1)'))
    complexIm(X1.*exp(1i*-angle(u*v1')),0,2,16));colorbar;set(gca,'ytick',[1 64],'xtick',1:4,'fontsize',16);...
    colormap hsv; freezeColors;title('Demodulated Theta');xlabel('Time (s)');

figure;subplot(4,1,2);imagesc(complexIm(v2,0,1));colorbar;axis off;
c = ceil(mod(angle(v2),2*pi)/(2*pi)*64);%(angle(v2*1i)+pi)
subplot(4,1,3);scatter(1:numel(v2),real(v2),[],angCol(c,:),'filled');
set(gca,'color','w','xtick',[],'ytick',[]);axis tight; colorbar;
%% panel 6
offSet = 1;
Xf1 = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).');...
[zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xf1 = [real(Xf1);imag(Xf1)];
Xf1 = zscore(Xf1,0,2);
t = conj(complex(W(p2,1:end/2),W(p2,end/2+1:end)))*complex(Xf1(1:65,:),Xf1(66:end,:));
%%
t = W(p2,:)*zscore(Xf1,0,2);
Wp = pinv(W);Wp = Wp(:,p2);
%%
[sp cellInfo] = hipSpikes('ec014.468',32/1.25);
spf = morFilter(sp,8,1250/32);clear sp;
s1 =  [35    17    19     7    26    27    60     6    33    39    50    54    18    67    21    31 14    36    46    44    40     1    56    59  24];
runTriggerViewSp(pos,v,spf(cellInfo == 1,:),[50 1],.05,s1(r(1:25)));
%% correlate high freq ICA and spiking
[sp cellInfo] = hipSpikes('ec014.468',32/1.25);
[A,tes,td] = allHighFreq('ec014.468.h5',1:8,pos,1);
tdf = morFilter(td,8,1250/32);
spf = morFilter(sp,8,1250/32);clear sp;
[~,~,teSp] = runTriggerViewSp(pos,v,spf(cellInfo >= 1,:),[50 1],.05);
[~,~,tes1] = runTriggerViewC(pos,v,tdf,[50 1],.05,eye(size(tdf,1)));
for i = 1:size(tes1,1)
    tes1f(i,:,:) = imfilter(squeeze(tes1(i,:,:)),fspecial('gaussian',5,1));
end
for i = 1:85
    teSpf(i,:,:) = imfilter(squeeze(teSp(i,:,:)),fspecial('gaussian',5,1));
end
ccf = abs(corr(tes1f(:,:)',teSpf(:,:)'));
figure;ccfa = ccf1;
for i = 1:18
[fx(i),fy(i)] = find(ccfa == max(ccfa(:)));
ccfa(fx(i),:) = 0;
ccfa(:,fy(i)) = 0;
subplot(6,6,2*i-1);imagesc(complexIm(imfilter(squeeze(tes1(fx(i),:,:)),fspecial('gaussian',5,1)),0,1));axis off;
subplot(6,6,2*i);imagesc(complexIm(imfilter(squeeze(teSp(fy(i),:,:)),fspecial('gaussian',5,1)),0,1));axis off;
end
%% compare different layers in 512-electrode data
%ss = [];
%for i =1:32
%x = getData('AB3-58.h5',(i-1)*16+(1:16),1000000);
%x = filtHigh(x,1250,100,8);
%s = std(x,0,2);
%ss = [ss s];
%end
%ss1 = ss(:);
%probes1= probes(:,[1:12 14 13 16 15]);
%for i = 1:size(probes,1)
%for j = 1:size(probes,2)
%ssa(i,j) = ss1(probes1(i,j)+1);
%end
%end
r = rippledet('AB3-59',512,154);
rm = mean(r,2);
rm = showGrid(rm,probes1);
rm = rm > graythresh(rm);
r = roipoly;
r1 = r & showGrid(1-rm,probes1);
%f = probes1(ssa < .11 & r)+1;
f = probes1(r1)+1;
fTest = zeros(512,1);
fTest(f) = 1;
figure;showGrid(fTest,probes1);
[A,W] = runTriggerICA(pos,v,Xf(f,:),.05);
%r1 = roipoly;
%f1 = probes1(ssa > .12 & r1)+1;
%[A1,W1] = runTriggerICA(pos,v,Xf(f1,:),.05);
%[A2,W2] = runTriggerICA(pos,v,Xf([f; f1],:),.05);
%% ica on anterior, posterior, and both shanks (AB3-58AP.mat)
[A2,W2,Z2] = runTriggerICA(pos,v,Xf(:,:),.05);
[A,W,Z] = runTriggerICA(pos,v,Xf(probes1(:) > 255,:),.05);
[A1,W1,Z1] = runTriggerICA(pos,v,Xf(probes1(:) <= 255,:),.05);
runTriggerViewC(pos,v,Xf(probes1(:) <= 255,:),[50 1],.05,W1);
runTriggerViewC(pos,v,Xf(probes1(:) > 255,:),[50 1],.05,W);
%% high freq ica in 2d environs
ccfa = corr(sp(cellInfo >= 1,:)',td');
for i = 1:18
[fx(i),fy(i)] = find(ccfa == max(ccfa(:)));
ccfa(fx(i),:) = 0;
ccfa(:,fy(i)) = 0;
end
td(1:2:35,:) = td(fy,:);
sp = sp(cellInfo >= 1,:);
td(2:2:36,:) = sp(fx,:);
A1(:,2:2:36) = A(:,fy);
A1(:,1:2:35) = A(:,fy);
for i = 1:16
subplot(6,6,2*i);title(sh(fx(i)));
end
%% cross frequency
[a b c d] = runTriggerViewBoth(X,Xf,v,pos,[.05 1],50,whiteningMatrix,phi,W);
 c2 = c;for i = 1:size(c,1)
c2(i,:,:) = imfilter(squeeze(c(i,:,:))./max(d,1),fspecial('gaussian',5,1));
c2(i,:,:) = c2(i,:,:)/max(c2(i,:));
end
a2 = a;for i = 1:size(a,1)
a2(i,:,:) = imfilter(squeeze(a(i,:,:))./max(b,1),fspecial('gaussian',5,1));
a2(i,:,:) = a2(i,:,:)/max(abs(a2(i,:)));
end
f1 = figure;f2 = figure;
for i = 1:64
figure(f1);subplot(8,8,i);plot((dewhiteningMatrix*squeeze(phi(:,i,:)))');
figure(f2);subplot(8,8,i);imagesc(squeeze(c2(i,:,:)));
end
cc = corr(abs(a2(:,:))',abs(c2(:,:))');
params.Fs = 1250/8;params.tapers = [1 1];params.pad = 2;
f1 = figure;f2 = figure;f3 = figure;cc1 = cc;for i = 1:20
[fx fy] = find(cc1 == max(cc1(:)));
cc1(fx,:) = 0;cc1(:,fy) = 0;
figure(f1);subplot(4,5,i);imagesc(abs([squeeze(c2(fy,:,:))/max(squeeze(c2(fy,:))); squeeze(a2(fx,:,:))/max(squeeze(a2(fx,:)))]));
figure(f2);subplot(4,5,i);set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(Z(:,fx),0,1)));
phiT = (dewhiteningMatrix*squeeze(phi(:,fy,:)));
plot(phiT');axis tight;
figure(f3);subplot(4,5,i);set(gca,'nextPlot','add','ColorOrder',squeeze(complexIm(Z(:,fx),0,1)));
[S,f] = mtspectrumc(phiT',params);
plot(f,S);axis tight;
end
%% sort positional indices
c = squeeze(mean(c,2));
posInds = find(max(abs(c),[],2) > 2);
[~,m] = max(abs(c(posInds,:)),[],2);
[~,s] = sort(m,'ascend');
posInds = posInds(s);
%% cross frequency 1
 u = u(:,1);
uc = angle(u);
uc = uc - min(uc);
uc = uc / (max(uc)+eps);
uc = 64-floor(uc*64);
cm = colormap;
cc = corr(abs(a2(:,:))',abs(c2(posInds,:))');
f1 = figure;f2 = figure;f3 = figure;
for i = 1:numel(posInds)
    figure(f1);subplot(4,5,i);set(gca,'nextPlot','add','ColorOrder',cm(uc,:));
    plot((dewhiteningMatrix*squeeze(phi(:,posInds(i),:)))');axis tight;
    figure(f2);subplot(4,5,i);set(gca,'nextPlot','add','ColorOrder',cm(uc,:));
    [~,m] = max(cc(:,i));
    plot(-(real(A(:,m)*exp(1i*linspace(0,2*pi,size(phi,3)))))');title(m);axis tight;
    figure(f3);subplot(4,5,i);imagesc(abs([squeeze(c2(posInds(i),:,:))/max(squeeze(c2(posInds(i),:)));...
        squeeze(a2(m,:,:))/max(squeeze(a2(m,:)))]));
end
%%4 examples for fig 4
ps = [9 10 16 18];
figure;
for i = 1:4
    subplot(2,4,i);set(gca,'nextPlot','add','ColorOrder',cm(uc,:));
    plot((dewhiteningMatrix*squeeze(phi(:,posInds(ps(i)),:)))');axis tight;
    [~,m] = max(cc(:,ps(i)));
    subplot(2,4,i+4);set(gca,'nextPlot','add','ColorOrder',cm(uc,:));
    plot(-(real(A(:,m)*exp(1i*linspace(0,2*pi,size(phi,3)))))');axis tight;
end
%%
[magAll timeAll numAll numSamps] = runConvAct(X,pos);
j = 1;
sPlot(squeeze(complex(numAll(j,:,:)/100,timeAll(j,:,:)./(max(1,numAll(j,:,:))))))
tempMag = squeeze(magAll(j,:,:));
[~,mx] = max(abs(tempMag),[],2);
tempTime = squeeze(timeAll(j,:,:)./max(1,numAll(j,:,:)));
for i = 1:size(tempMag,1)
tM(i,:) = circshift(tempMag(i,:),[0 50-mx(i)]);
tT(i,:) = circshift(tempTime(i,:),[0 50-mx(i)]);
end
plot(tT')
%%FIG 1 anatomy
figure;imagesc(pos(:,:,3));colormap gray;axis image;
[xs,ys] = meshgrid((1:32)*50,(1:8)*300);
scale = 5/3;
hold all;scatter(ys(:)*scale-150,xs(:)*scale,'r','filled')
hold all;scatter(ys(:)*scale-150+4000,xs(:)*scale,'r','filled')
%load elecPos.mat
x8 = reshape(cc1(1,:),[8 8]);
y8 = reshape(cc1(2,:),[8 8]);
y8 = y8/mean(mean(diff(y8)))*20;
x8 = x8/mean(mean(diff(x8')))*200;
scatter(x8(:)*scale+1200,-y8(:)*scale,'y','filled');
x8 = x8(:,5:8);y8 = y8(:,5:8);
scatter(x8(:)*scale+1200,-y8(:)*scale,'g','filled');
scatter(x8(:)*scale+1200,-y8(:)*scale,'y');
%% find sequences of ic activations
allCorr = zeros(numel(lags),512,512);for j = 1:numel(lags)
allCorr(j,:,:) = Xfd*circshift(Xfd,[0 lags(j)])';
end
[mVal,mInds] = max(abs(allCorr));
mVal = squeeze(mVal);mInds = squeeze(mInds);

%%try 2
for i = 1:512
m(i) = max(bwlabel(bwmorph(abs(Xfd(i,:)) > 4,'dilate',10)));
end
figure;plot(m)
fi = find(m > 10);
%sPlot(abs(Xfd(fi,:)))
%[inds,fi] = sort(m,'descend');
xy = [];for i = 1:numel(fi)
temp = imfilter(squeeze(actMax(fi(i),:,:)),fspecial('gaussian',5,1));
[xy(i,1) xy(i,2)] = find(temp == max(temp(:)));
end
[a b c] = mtspo_ga(xy,1-mVal(fi,fi),1);
%%plot vel plus sequences of activation in maze
vel = angVel(pos);
veld = [];
for i = 1:2
veld(:,i) = decimate(vel(:,i),8);
end
figure;subplot('position',[0 .8 1 .2]);
inds = 5900:6500;
veld = filtLow(veld(inds,1),1250/32/8,1);
plot(veld,'b','linewidth',2);hold all;
vSig = filtLow(v1(inds),1250/32/8,1);
vSig = vSig - prctile(vSig,5);
b = veld > 2;
for i = 1:4
    b = b - bwmorph(b,'endpoints');
end
b = bwmorph(b,'dilate',7);
for i = 1:3
    b = b - bwmorph(b,'endpoints');
end
hold all;
plot(vSig/max(vSig)*8,'r','linewidth',2);axis tight;
set(gca,'xtick',[],'ytick',[]);
subplot('position',[0 0 1 .8]);hold all;
imagesc(filtLow(abs(Xfd(fi,inds)),1250/32/8,1));
b = find(bwmorph(b,'endpoints'));
for i = 1:numel(b)
    subplot('position',[0 0 1 .8]);plot(b(i)*ones(2,1),[0 numel(fi)],'w--','linewidth',2);
    %subplot('position',[0 .8 1 .2]);plot(b(i)*ones(2,1),[-.5 10],'k--','linewidth',2);
end
subplot('position',[0 .8 1 .2]);set(gca,'ylim',[-.5 max(veld)*1.3]);
%% make patchwork for maze activations
[im frameCol] = superImp(absMean,fi,1,5,6);
%% zoom into home base
[a b] = restlessBox(Xfd(fi,:),pos,[-inf 160 160 300],frameCol,37:39,33:36);
%% draw boundaries on array
for i = 1:3
xx = min(bounds1{i}.xs):max(bounds{i}.xs);
yy = interp1(bounds1{i}.xs,bounds{i}.ys,xx,'cubic');
plot(xx,yy,'w--','linewidth',2);
end