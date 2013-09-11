function phaseFlow(phi,ind)

sz = [8 8];
[xs ys] = meshgrid(1:sz(1),1:sz(2));
grid = [xs(:); ys(:)] + .5;
per = 8/1250;
phiT = squeeze(phi(:,ind,:));
vals = 1./(per*(2:2:40));
wav = 'cmor5-1';
f = scal2frq(vals,wav,per);
t = zeros(size(phiT,1),numel(vals),size(phiT,2));
for i = 1:size(phiT,1)
    t(i,:,:) = cwt(phiT(i,:),vals,wav);
%    subplot(8,8,i);imagesc(complexIm(t(i,:,:),0,1));
end
ts = squeeze(mean(abs(t),3));
figure;subplot(211);plot(ts');
params.Fs = 1/per; params.tapers = [1 1]; params.pad = 2;
[S,f1] = mtspectrumc(phiT',params);
subplot(212);plot(f1,sqrt(S));
ts = mean(ts);
c= colormap(jet(size(phi,3)));
xdim = ceil(sqrt(numel(f)));ydim = ceil(numel(f)/xdim);
f1 = figure;
f2 = figure;f3 = figure;
scale = 5;
for i = 1:numel(vals)
    tt = squeeze(t(:,i,:));
    [u,~,v] = svds(tt,1);
    va = angle(mean(conj(u))*v);
    tt = bsxfun(@times,tt,exp(1i*va).');
    figure(f1);subplot(xdim,ydim,i);imagesc(complexIm(tt,0,1));
    figure(f2);subplot(xdim,ydim,i);imagesc(complexIm(flipud(reshape(u*exp(1i*angle(mean(conj(u)))),[8 8])),0,1,60/f(i)));
    figure(f3);subplot(xdim,ydim,i);
    for j = 1:size(phi,3)
        [dx dy] = angGradient(reshape(tt(:,j),sz));
        dx = dx/f(i)*max(f);
        dy = dy/f(i)*max(f);
        hold all;
%        scatter(scale(1)*dx(:)+xs(:),scale(end)*dy(:) + ys(:),5,c(j,:),'filled');
        quiver(dx,dy,'color',c(j,:));
    end
    title([f(i) ts(i)]);
    axis tight;
%    set(gca,'xlim',[0 9],'ylim',[0 9]);
end