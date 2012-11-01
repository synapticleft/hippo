function phaseEvol(dat,u,s,v,inds,mName,sp,c)

trails = 0;lags = 10;
rdim = 2;
skip = 1;
if ~exist('inds','var') || isempty(inds)
    inds = 2:skip:size(dat,2);
end
data1 = u(:,1:rdim)*s(1:rdim,1:rdim)*v(inds,1:rdim)';
[~,ord] = sort(circ_dist(angle(u(:,1)),circ_mean(angle(u(:,1)))),'descend');
%u(:,1) = 1;

h = figure;

if ~exist('c','var')
    c = getCol(size(dat,1));
else
    c = c(ord);
end
if exist('sp','var') && ~isempty(sp)
    sp(sp == -1) = nan;
    if numel(sp) < size(v,1)
        v((numel(sp)+1):end,:) = [];
        dat(:,(numel(sp)+1):end) = [];
    else
        sp((size(v,1)+1):end) = [];
    end
    if any(isnan(sp))
        nanInds = find(~isnan(sp));
        sp = interp1(nanInds,sp(nanInds),1:numel(sp));
        nanInds = isnan(sp');
        sp(nanInds) = [];
        v(nanInds,:) = [];
        dat(:,nanInds) = [];
        sp = sp - mean(sp);
        sp = sp/std(sp);
    end
    %s = std(sp(:));
    c1 = getCol(size(sp,1));
    [~,ind] = sort(mean(abs(sp),2),'ascend');
    c1(ind,:) = c1;
else
    c1 = c;
end
circBuff1 = nan*ones(size(c1,1),lags);circBuff2 = nan*ones(size(dat,1),lags);
%inst = circ_mean(circ_dist(angle(v(1:(end-1),1)),angle(v(2:end,1))));
md = -circ_mean(circ_dist(angle(v(2:end,1)),angle(v(1:(end-1),1))));%circ_diff(v(:,1).')');
for i = inds
    figure(h);
phaseShift = 0;%-angle(mean(u(:,1))*v(i-1,1)');%-i*md*skip;%-angle(u(:,1)*v(i,1)');%angle(u(:,1)*circ_mean(angle(dat(:,i)),abs(dat(:,i))));%
if exist('sp','var') && ~isempty(sp)
    %imagesc([-.1 .1],[-2 2],[real(sp(:,i)); spa(i-inds(1)+1)] ,s*[-1 1]*2);%colormap gray;
    ps = phaseShift;%-angle(v(i,1)')
    temp = sp(:,i).*exp(1i*ps);
else
    temp = dat(ord,i).*exp(1i*phaseShift);
%    temp = temp/mean(abs(temp));
end
scatter(real(temp),imag(temp),60,'k','p','filled');axis image;%c1
hold on;
%scatter(real(temp(reds)),imag(temp(reds)),60,'r','p','filled');
circBuff1 = circshift(circBuff1,[0 1]);
circBuff1(:,1) = temp;
hold on;
temp = data1(ord,i-inds(1)+1).*exp(1i*phaseShift);
%%%
scatter(real(temp),imag(temp),30,c,'o','filled');%c
%tCol(ord,:) = c;
%figure;image(reshape(tCol,[16 6 3] ));
circBuff2 = circshift(circBuff2,[0 1]);
circBuff2(:,1) = temp;
if trails
    set(gca,'ColorOrder',c1);plot(circBuff1.');
    set(gca,'ColorOrder',c);plot(circBuff2.');
    %scatter(real(circBuff1(:)), imag(circBuff1(:)),c1);
    %scatter(real(circBuff2(:)), imag(circBuff2(:)),c);
end
temp = mean(abs(dat(:,i)))*exp(1i*(-3.2:.1:3.2));
plot(real(temp),imag(temp),'k');
%temp = mean(abs(dat(:,i)))*exp(1i*-(angle(v(i,1))))*exp(1i*phaseShift);%+inst*i
temp = s(1)*abs(mean(u(:,1)))*conj(v(i,1))*exp(1i*phaseShift);%exp(1i*-(angle(v(i,1))))
%scatter(real(temp),imag(temp),50,'k','filled');
temp = 10*s(1)*abs(mean(u(:,2)))*conj(v(i,2))*exp(1i*phaseShift);%exp(1i*-(angle(v(i,1))))
%scatter(real(temp),imag(temp),50,'r','filled');
hold off;
set(gca,'xlim',[-2 2],'ylim',[-2 2],'xticklabels',[],'yticklabels',[],'fontsize',16);
title(num2str((i-inds(1)+1)*32/1250));
m(i-inds(1)+1) = getframe(gcf);
drawnow;
end
if exist('mName','var')
    movie2avi(m,mName,'fps',8);
end
%c = colormap;
%mpgwrite(m,c,'spin.avi');


function c = getCol(n)
c = repmat(linspace(0,1,n)',[1 3]);
for i = 1:3
c(:,i) = max(0,min((1-abs(c(:,i)-i/4)*2.5)*1.5,1));
end
c = fliplr(c);