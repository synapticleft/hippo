function posVsTimeScatter(pos,v,bounds,inds)%,sp

if size(v,2) > size(v,1)
    v = v.';
end
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),1);
end
nanInds = find(~isnan(pos));
pos = interp1(nanInds,pos(nanInds),1:numel(pos));
nanInds = isnan(pos);
pos = pos(~nanInds);v = v(~nanInds,:);
if exist('inds','var')
    pos = pos(end-inds:end);%(1:inds);
    v = v(end-inds:end,:);%(1:inds,:);
end
b = nan*ones(numel(pos),1);
b(pos < bounds(1)) = -1;b(pos > bounds(2)) = 1;
figure;plot(b);
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:numel(pos));
b = [0 diff(b)];
angCol = colormap('hsv');absCol = colormap('jet');
v = bsxfun(@rdivide,v,std(v));
vp11 = v(2:end,1).*conj(v(1:end-1,1))./abs(v(1:end-1,1));
vp11 = [0; vp11];
vp12 = v(:,2);%v(:,1).*conj(v(:,2))./abs(v(:,1));%
%vp12 = log(vp12);
% state = [real(vp12) imag(vp12)];
% [data al] = binData(state(b>0,:),pos(b>0),[20 20]);
% figure;subplot(221);imagesc(data);
% subplot(222);imagesc(al);
% [data al] = binData(state(b<0,:),pos(b<0),[20 20]);
% subplot(223);imagesc(data);
% subplot(224);imagesc(al);
% [data al] = binData(pos(b<0).',vp12(b<0).',[20]);
% figure;plot(real(data));hold all;plot(imag(data));plot(al);
% return
figure;subplot(211);scatter(1:sum(b>0),pos(b>0),50*abs(vp12(b>0)),angCol(phase2Col(angle(vp12(b>0))),:),'filled');axis tight;
subplot(212);scatter(1:sum(b<0),pos(b<0),50*abs(vp12(b<0)),angCol(phase2Col(angle(vp12(b<0))),:),'filled');axis tight;
%vp12 = vp12 + 1i*(1:numel(vp12))'/200;%vp12 = log(vp12);
pos = pos - min(pos);pos = pos / max(pos);pos = 1-pos;
figure;subplot(121);%hold all;plot(real(vp12(b<0)),imag(vp12(b<0)));
scatter(real(vp12(b<0)),imag(vp12(b<0)),15,absCol(abs2Col(pos(b<0)),:),'filled');
pos = 1-pos;
subplot(122);%hold all;plot(real(vp12(b>0)),imag(vp12(b>0)));
scatter(real(vp12(b>0)),imag(vp12(b>0)),15,absCol(abs2Col(pos(b>0)),:),'filled');
drawnow;

function c = phase2Col(ang)
c = ceil((ang+pi)/(2*pi)*64);


function c = abs2Col(a)
if isnan(a)
    a = 0;
end
c = 1 + floor(64*min(.99,a));

function [data al] = binData(x,y,nbins)
[size(x) size(y)]
for i = 1:size(x,2)
    s = std(x(:,i));m = mean(x(:,i));
    x(:,i) = min(max(x(:,i),m-2*s),m+2*s);
    x(:,i) = x(:,i) - min(x(:,i));
    x(:,i) = ceil(nbins(min(i,numel(nbins)))*max(eps,x(:,i))/max(x(:,i)));
    x(:,i) = min(nbins(min(i,numel(nbins))),x(:,i));
end
figure;plot(x);
data = accumarray(x,y,nbins,@mean);
al = accumarray(x,y,nbins,@std);
data(data == 0) =nan;