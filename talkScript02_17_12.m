%%1st PC props sm96
temp = (reshape(angle(u(:,1)),[16 6]));
figure;subplot(211);
imagesc(temp);
[xt yt] = myGradient(temp);xt = -xt;yt = -yt;
hold all;quiver(xt,yt,'w','linewidth',2);
set(gca,'fontsize',16)
title 'angle 1st PC'
subplot(212);
drawnow;
X = getData('sm9603Circ1.h5',10000,0,1);%'ec013.670.h5'
params.Fs = 1250;params.tapers = [3 5];
%for i= 1:96
%[S,f] = mtspectrumc(X(i,:),params);
%f = find(f > 200 & f < 500);
%pyr(i) = mean(S(f));
%end
imagesc((reshape(pyr,[16 6])));
hold all;quiver(xt,yt,'w','linewidth',2);
set(gca,'fontsize',16)
title 'power 200-500 Hz';
%% 1st PC props ec014
pyr = [];
temp = flipud(reshape(angle(u(:,1)),[8 8]));
figure;subplot(211);
imagesc(temp);
[xt yt] = myGradient(temp);xt = -xt;yt = -yt;
hold all;quiver(xt,yt,'w','linewidth',2);
set(gca,'fontsize',16)
title 'angle 1st PC';
colorbar
subplot(212);
drawnow;
X = getData('ec014.277.h5',10000,0,1);%'ec013.670.h5'
params.Fs = 1250;params.tapers = [3 5];
for i= 1:64
[S,f] = mtspectrumc(X(i,:),params);
f = find(f > 200 & f < 500);
pyr(i) = mean(S(f));
end
imagesc(flipud(reshape(pyr,[8 8])));
hold all;quiver(xt,yt,'w','linewidth',2);
set(gca,'fontsize',16)
title 'power 200-500 Hz';colorbar
%% 2nd fig
figure;subplot(221);imagesc(angle(reshape(u(:,1),[8 8])));
set(gca,'fontsize',16)
title 'angle 1st PC'
subplot(223);imagesc(real(reshape(u(:,2),[8 8])))
set(gca,'fontsize',16)
title 'Re(2nd PC)'
temp0 = u(:,1)*s(1,1)*v(:,1)';
for i = 1:64
temp(i) = corr(Xf(i,:).',temp0(i,:).');
end
subplot(222);imagesc(reshape(abs(temp),[8 8]).^2,[.83 1]);
title 'variance explained, 1PC'
colorbar
for i = 1:64
temp(i) = corr(Xf(i,:).',temp0(i,:).');
end
 subplot(224);imagesc(reshape(abs(temp),[8 8]).^2,[.83 1]);
set(gca,'fontsize',16)
title 'variance explained, 2PCs'
colorbar
%% 3rd
inds = 191500+(1:205);%86000:90000;
temp = getData('ec014.277.h5',205*32,191500*32);
figure;subplot(511);plot((1:numel(inds))*32/1250,vel(inds,1));axis tight;%sqrt(vel(inds,1).^2)+sqrt(vel(inds,2).^2));
set(gca,'fontsize',16,'xtick',[]);ylabel 'velocity';colorbar
subplot(512);imagesc((1:numel(inds)*32)/1250,1:64,temp,[prctile(temp(:),1) prctile(temp(:),99)]);colorbar%real(u(:,1)*s(1,1)*v(inds,1)'));colorbar
set(gca,'fontsize',16);axis off;%title 'LFP';
%subplot(413);imagesc(abs(u(:,1)*s(1,1)*v(inds,1)'));
temp = imag(bsxfun(@times,Xf(:,inds),conj(v(inds,1)./abs(v(inds,1)))'));
subplot(513);imagesc((1:numel(inds))*32/1250,1:64,temp,[prctile(temp(:),1) prctile(temp(:),99)]);colorbar
set(gca,'fontsize',16);axis off;%title 'Im(projected LFP)'
subplot(514);imagesc((1:numel(inds))*32/1250,1:64,temp-imag(bsxfun(@times,u(:,1)*s(1,1)*v(inds,1)',conj(v(inds,1)./abs(v(inds,1)))')),[-1 1]/2);colorbar
set(gca,'fontsize',16);axis off;%title 'Im(2-PC approximation)'
subplot(515);imagesc((1:numel(inds))*32/1250,1:64,temp-imag(bsxfun(@times,u(:,1:2)*s(1:2,1:2)*v(inds,1:2)',conj(v(inds,1)./abs(v(inds,1)))')),[-1 1]/2);colorbar
set(gca,'fontsize',16,'ytick',[]);%title 'Im(2-PC approximation)'
%% 4th
runTrigger(pos,v(:,1:2),.3,200);
%% 5th
for i = 1:4
ccs16(i,:) = posVsSpike(pos,spf(inds(1:thisMany(i)),:),v,[.1 .9],[60]);
end
figure;plot(ccs13(:,[1 3]),'b--','linewidth',2);
hold all;plot(ccs13(:,[2 4]),'b','linewidth',2);
hold all;plot(ccs14(:,[2 4]),'r','linewidth',2);
hold all;plot(ccs14(:,[1 3]),'r--','linewidth',2);
hold all;plot(ccs16(:,[1 3]),'k--','linewidth',2);
hold all;plot(ccs16(:,[2 4]),'k','linewidth',2);
set(gca,'fontsize',16)
set(gca,'xticklabel',[5 10 20 40]);
set(gca,'xtick',[1 2 3 4],'xticklabel',[5 10 20 40]);
xlabel '# neurons'
ylabel 'R squared'
title 'PC2 prediction by neuron firing'