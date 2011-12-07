X = getData(20000,5);

%%for spsvd
[u,s,v] = svds(X,20);
v = s*v';
params.tapers = [7 13];
params.Fs = 250;
params.fpass = [0 70];
[sv,sp,fm] = spsvd(v(:,1:3000)',params);
sPlot(sv');
for i = 1:inf
inds = 0:2:floor(max(f));inds = inds/ceil(max(f))*numel(f);inds = 1+floor(inds);
for j = 1:4
subplot(2,2,j);showGrid(real(bsxfun(@times,u*sp(:,inds,j),exp(1i*i/20))),[8 8],[1 1]);
end
end

%%For SSA
[a,b,c,d] = ssa1(X,50,20);
SSAMovies(d,b,[8 8],[1 1]);
[S,f] = mtspectrumc(c,params);
ch = hilbert(c.');
for i = 1:10
for j = 1:10
    subplot(10,10,10*(i-1)+j);imagesc(hist3([angle(ch(i,:));angle(ch(j,:))]',[100 100]));axis off;
end
end
[u1 s1 v1] = svds(ch,50);
[Sh,fh] = mtspectrumc(v1,params);
sPlot(Sh',fh);
v1 = v1.';for i = 1:10
for j = 1:10
    subplot(10,10,10*(i-1)+j);imagesc(hist3([angle(v1(i,:));angle(v1(j,:))]',[100 100]));axis off;
end
end

%for AR
[u,s,v] = svds(X,10);
v = s*v';
[w,a,c] = arfit(v',1,5,'zero');
Xar1 = arsim(w,a,c,10000)';
Xar1 = u*Xar1;
[Sar,fa] = mtspectrumc(Xar1(1,:),params);
[S,f] = mtspectrumc(X(1,:),params);
figure;plot(f,S);hold all;plot(fa,Sar);

Xf = hipfilter(X,5,11,250,8);
Xar1f = hipfilter(Xar1,5,11,250,8);
figure;for i = 1:inf
    subplot(221);imagesc(reshape(X(:,i),[8 8]));
    subplot(222);imagesc(reshape(Xar1(:,i),[8 8]));
    subplot(223);imagesc(reshape(Xf(:,i),[8 8]));
    subplot(224);imagesc(reshape(Xar1f(:,i),[8 8]));
    pause(.05);
    drawnow;
end

[psi,x] = cmorwavf(-1,1,100,1/500,8.4);
for i = 1:64
Xf(i,:) = fliplr(filter(conj(psi),1,fliplr(filter(psi,1,X(i,:)))));
end
Xm = circ_mean(angle(Xf),abs(Xf));
Xh1 = bsxfun(@times,Xh,exp(-1i*Xm));
figure;imagesc(reshape(circ_mean(angle(Xh1),[],2),[8 8]))
figure;imagesc(reshape(mean(abs(Xf(:,(3050*50):(3400*50))),2),[8 8]))
figure;imagesc(reshape(circ_mean(angle(Xh1(:,(3050*50):(3400*50))),[],2),[8 8]),[-.25 .5])
figure;imagesc(reshape(circ_std(angle(Xh1),[],[],2),[8 8]))
mean(circ_dist(angle(Xf(1,(3050*50):(3400*50))),angle(Xf(1,((3050*50):(50*3400))-1))))
1900-2350
%angular corr
[u s v] = svds(Xf(:,1:100000),1);
Xi = u*s*v';
[u s v] = svds(Xf(:,1:100000),2);
Xii = u*s*v';
[u s v] = svds(Xf(:,1:100000),3);
Xiii = u*s*v';
for i =1:64
c(1,i) = circ_corrcc(angle(Xi(i,:)),angle(Xf(i,1:100000)));
c(2,i) = circ_corrcc(angle(Xii(i,:)),angle(Xf(i,1:100000)));
c(3,i) = circ_corrcc(angle(Xiii(i,:)),angle(Xf(i,1:100000)));
end
figure;for i = 1:3
subplot(3,1,i);imagesc(reshape(c(i,:),[8 8]),[.9 1]);
end
figure;imagesc(circ_dist(angle(Xf(:,1000:2000)),angle(Xf(:,1000:2000)-Xi(:,1000:2000))),[-pi pi]);colormap hsv
temp = circ_dist(angle(v(:,1)),angle(v(:,2)));
mean(abs(circ_dist(temp(1:(end-1)),temp(2:end))))
mean(circ_dist(temp(1:(end-1)),temp(2:end)))
%how much of angular spread does 2nd component explain
sf = circ_std(angle(Xf),abs(Xf));
[u s v] = svds(Xf(:,1:100000),2);
sf2 = circ_std(angle(u*s*v'),abs(u*s*v'));
corr(sf',sf2')
%%%%%%%%%%%%%%%DISPLAY GRADIENT
figure;for i = 1:inf
scatter(circ_dist(angle(u(:,1)),circ_mean(angle(u(:,1))))*abs(v(i,1)),circ_dist(angle(Xf(:,i)),angle(u(:,1)*v(i,1)'))*abs(v(i,1)));hold on;
scatter(circ_dist(angle(u(:,1)),circ_mean(angle(u(:,1))))*abs(v(i,1)),circ_dist(angle(Xf1(:,i)),angle(u(:,1)*v(i,1)'))*abs(v(i,1)),'r');
hold off;
set(gca,'xlim',[-1 2]/1000,'ylim',[-3 3]/1000);drawnow;
end
%%%%%%%%scatter phases
figure;lag = 0;for i = 1000:inf
subplot(121);
temp = data1(:,i+lag).*exp(1i*angle(v(i,1)));
scatter(real(temp),imag(temp),20,c,'o','filled');hold on;
temp = mean(abs(Xf(:,i+lag)))*exp(1i*(-pi:.1:pi));plot(real(temp),imag(temp),'k');
hold off;
set(gca,'xlim',[-2 2],'ylim',[-2 2]);
subplot(122);
temp = data1(:,i+lag);
scatter(real(temp),imag(temp),20,c,'o','filled');hold on;
temp = mean(abs(Xf(:,i+lag)))*exp(1i*(-pi:.1:pi));plot(real(temp),imag(temp),'k');
hold off;
set(gca,'xlim',[-2 2],'ylim',[-2 2]);drawnow;
end
%%OR
posa(posa == -1) = nan;
for i = 1:4
posb(:,i) = resample2((1:10000)/50,posa(:,i),(1:10000)/39.06);
end
subplot(122);
scatter(posb(i+lag1,1),posb(i+lag1,2),'b','filled');
hold on;scatter(posb(i+lag1,3),posb(i+lag1,4),'r','filled');
plot(posb(i:(i+lag1),1),posb(i:(i+lag1),2),'b');
plot(posb(i:(i+lag1),3),posb(i:(i+lag1),4),'r');
set(gca,'xlim',[0 250],'ylim',[0 250]);
hold off;drawnow;
%%%Spike triggered histograms
bins{1} = -pi:.05:pi;bins{2} = bins{1};
im1 = hist3([angle(v(:,1)) circ_dist(angle(v(:,1),angle(v(:,2))))],bins);
 figure;for i = 1:36
subplot(6,6,i);
x = v(find(sp(inds(i),:)),:);
imagesc(hist3([angle(x(:,1)) circ_dist(angle(x(:,1)),angle(x(:,2)))],bins)./im1)
axis image off;
end
%%%%%%%%%%%