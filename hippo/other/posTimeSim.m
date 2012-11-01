function posTimeSim(pos,v,bounds,gauss)
shift = 1;
accumbins = gauss(1);pad = 0;
%[gx gy] = meshgrid(1:gauss(1));% = (1:accumbins)*(1:accumbins)';
%gaussWeights = exp(-((gx-gauss(2)).^2+(gy-gauss(2)).^2)/gauss(3).^2);
%gaussWeights = exp(-((1:gauss(1))-gauss(2)).^2/gauss(3).^2);
if size(v,2) > size(v,1)
    v = v.';
end
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
nanInds = find(~isnan(pos(:,1)));
pos(:,1) = interp1(nanInds,pos(nanInds,1),1:size(pos,1));
pos(:,2) = interp1(nanInds,pos(nanInds,2),1:size(pos,1));
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);
 pos = bsxfun(@minus,pos,mean(pos));%pos = bsxfun(@rdivide,pos,std(pos));
 [a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
v = bsxfun(@rdivide,v,std(v));
vp11 = v((shift+1):end,1).*conj(v(1:end-shift,1))./abs(v(1:end-shift,1));
vp11 = [zeros(shift,1); vp11];
vp12 = v(:,1).*conj(v(:,2))./abs(v(:,1));vp12 = vp12 - mean(vp12);vp12 = complex(real(vp12)/std(real(vp12)),imag(vp12)/std(imag(vp12)));
%vp12 = vp11;
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
lrRuns = bwlabel(b>0);rlRuns = bwlabel(b<0);
lrAvg = zeros(accumbins+2*pad,max(lrRuns));
rlAvg = zeros(accumbins+2*pad,max(rlRuns));
medBinTime = round((median(hist(rlRuns,1:max(rlRuns)))+median(hist(lrRuns,1:max(lrRuns))))/(2*accumbins));
lrTimes = zeros(accumbins,max(lrRuns));
rlTimes = zeros(accumbins,max(rlRuns));
binsTimes = bounds(1)+((1:accumbins)-.5)/accumbins*diff(bounds);
accumbins = accumbins+pad*2;
bins = (bounds(1)-pad)+((1:accumbins)-.5)/accumbins*(diff(bounds)+2*pad);
for i = 1:max(lrRuns)
    inds = find(lrRuns == i);inds = (min(inds)-pad*medBinTime):max(inds)+pad*medBinTime;
    lrAvg(:,i) = csaps(pos(inds,1),vp12(inds),1-1e-7,bins);
%    lrTimes(:,i) = csaps(pos(inds,1),inds-min(inds),1-1e-9,binsTimes);%,pos(inds,1),1-1e-7
end
%drawnow;
%imagesc(abs(cov(lrAvg')));return
%axis tight;subplot(122);hold all;
for i = 1:max(rlRuns)
    inds = find(rlRuns == i);inds = (min(inds)-pad*medBinTime):max(inds)+pad*medBinTime;
    rlAvg(:,i) = csaps(pos(inds,1),vp12(inds),1-1e-7,bins);
%    rlTimes(:,i) = csaps(pos(inds,1),inds-min(inds),1-1e-9,binsTimes);%,pos(inds,1),1-1e-7
end
%lrAvg = bsxfun(@times,lrAvg,gaussWeights');
%rlAvg = bsxfun(@times,rlAvg,gaussWeights');
%[~,indsrl] = sort(mdscale(pdist([real(rlAvg); imag(rlAvg)]','correlation'),1));
figure;title('projection PC2->PC1','fontsize',16);
subplot(223);imagesc(angle(rlAvg));colormap hsv;%(:,indsrl)
[~,indslr] = sort(mdscale(pdist([real(lrAvg); imag(lrAvg)]','correlation'),1));
set(gca,'fontsize',16);
xlabel 'trial number'
ylabel 'position'
subplot(221);imagesc(angle(lrAvg));colormap hsv;%(:,indslr)
%figure;subplot(212);imagesc(abs(rlAvg));%(:,indsrl)
%subplot(211);imagesc(abs(lrAvg));%(:,indslr)
%figure;subplot(212);imagesc(real(rlAvg));%(:,indsrl)
%subplot(211);imagesc(real(lrAvg));%(:,indslr)
rlVa(:,:,1) = (angle(rlAvg)+pi)/(2*pi);rlVa(:,:,2) = 1;
rlVa(:,:,3) = power(abs(rlAvg)/max(abs(rlAvg(:))),.5);
rlVa = hsv2rgb(rlVa);
lrVa(:,:,1) = (angle(lrAvg)+pi)/(2*pi);lrVa(:,:,2) = 1;
lrVa(:,:,3) = power(abs(lrAvg)/max(abs(lrAvg(:))),.5);
lrVa = hsv2rgb(lrVa);
subplot(222);image(lrVa);subplot(224);image(rlVa);
set(gca,'fontsize',16);
xlabel 'trial number'
return
%axis tight;
subplot(222);plot(real(mean(lrAvg,2)));hold all;plot(-imag(mean(lrAvg,2)));
plot(abs(mean(lrAvg,2)));plot(-angle(mean(lrAvg,2)));%plot(std(lrAvg,0,2));axis tight;
subplot(224);plot(real(mean(rlAvg,2)));hold all;plot(-imag(mean(rlAvg,2)));
plot(abs(mean(rlAvg,2)));plot(-angle(mean(rlAvg,2)));%plot(std(rlAvg,0,2));axis tight;
%figure;plot(max(lrTimes));hold all;plot(max(rlTimes));
lrTimes = diff(lrTimes);rlTimes = diff(rlTimes);
%lrTimes = bsxfun(@minus,lrTimes,mean(lrTimes,2));rlTimes = bsxfun(@minus,rlTimes,mean(rlTimes,2));
figure;subplot(211);imagesc(angle(lrAvg));
subplot(212);imagesc(angle(rlAvg));colormap hsv;
figure;subplot(211);imagesc(abs(lrAvg));
subplot(212);imagesc(abs(rlAvg));
figure;subplot(211);imagesc(lrTimes,mean(lrTimes(:))+std(lrTimes(:))*[-1 1]*2);
subplot(212);imagesc(rlTimes,mean(rlTimes(:))+std(rlTimes(:))*[-1 1]*2);