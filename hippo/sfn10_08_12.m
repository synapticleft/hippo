%% Demodulation demo
nSteps = 1000;
x = linspace(0,20*pi,nSteps);
y = exp(1i*x);
A = [zeros(200,1); ones(300,1); zeros(500,1)]/2+1;
win = 20;
A = filtfilt(gausswin(win),sum(gausswin(win)),A);
dP = [zeros(600,1); ones(300,1); zeros(100,1)]*-pi/2;
dP = filtfilt(gausswin(win),sum(gausswin(win)),dP);
z = A'.*exp(1i*(x+dP'));zd = z.*conj(y);
angCol = colormap('hsv');
c = ceil((angle(z)+pi)/(2*pi)*64);
figure;plot3(1:nSteps,real(y),imag(y),'k');
hold all;scatter3(1:nSteps,real(z),imag(z),[],bsxfun(@times,angCol(c,:),(A/max(A))),'filled');
sh = 4;
plot3(1:nSteps,imag(y),real(y)-sh,'k');
c = ceil((angle(zd)+pi)/(2*pi)*64);
hold all;scatter3(1:nSteps,imag(zd),real(zd)-sh,[],bsxfun(@times,angCol(c,:),(A/max(A))),'filled');
set(gca,'xtick',[],'ytick',[],'ztick',[],'linewidth',2);
%%
inds = 52965:53140;
%temp = filtLow(angVel(pos)',1250/32,2);
[posa,s,u] = svds(pos(:,1:2),1);
posa = s*posa;
temp = filtLow(diff(posa(inds)),1250/32,2);
figure;subplot(411);plot(temp*1250/32,'k','linewidth',2);axis tight;
colorbar;set(gca,'xtick',[],'fontsize',16);title('Velocity');ylabel('cm/s');
indsa = inds(1)*4:inds(end)*4;
X1 = morFilter(X(:,indsa(1)-1000:indsa(end)+1000),8,1250/8);X1 = X1(:,1001:end-1000);
[u,~,v1] = svds(X1,1);
subplot(412);imagesc(X(:,indsa));set(gca,'ytick',[1 64],'xtick',[],'fontsize',16);colormap jet;colorbar;freezeColors;
subplot(413);imagesc(complexIm(X1,0,1));colorbar;set(gca,'ytick',[1 64],'xtick',[],'fontsize',16);
subplot(414);imagesc(linspace(0,size(X1,2)/1250*8,size(X1,2)),1:64,...bsxfun(@times,X1,exp(1i*angle(v1)'))
    complexIm(X1.*exp(1i*-angle(u*v1')),0,2,16));colorbar;set(gca,'ytick',[1 64],'fontsize',16);colormap hsv; freezeColors;