function stateHists2(pos,resp)
%2D histogram of neural response properties
%   Detailed explanation goes here

warning off all;
pos(pos == -1) = nan;
v = diff(pos);
%v = angVel(pos);
v = filtLow(v',1250/32,3)';
pos(1,:) = [];
%state = [pos(:,1) v(:,1)];
vp12 = resp(:,1).*conj(resp(:,2));vp12 = filtLow(vp12,1250/32,1);
vp1 = resp(:,1).*conj(resp(:,1));vp1 = filtLow(vp1,1250/32,1);
vp2 = resp(:,2).*conj(resp(:,2));vp2 = filtLow(vp2,1250/32,1);
vp11 = resp(1:end-1,1).*conj(resp(2:end,1));vp11 = filtLow(vp11,1250/32,1);
vp11 = [0;vp11];%vp11 = gsorth(vp11);
vp12 = vp12;%./vp1
state = [real(vp12) imag(vp12)];
%state = [real(vp11) imag(vp11)];
figure;plot(v(:,1)./max(v(:,1)));hold all;
%plot(v(:,2)./max(v(:,1)));
plot(bsxfun(@rdivide,state,max(eps,max(state(:)))));
%figure;imagesc(log(hist3(state,[100 100])));
[data al] = binData(state,v(:,1),[40 40]);%resp
figure;h = imagesc(data,[0 1.5]);
figure;imagesc(log(al));
drawnow;
%set(gca,'color','k');
%set(h,'AlphaData',al);

function [data al] = binData(x,y,nbins)
for i = 1:size(x,2)
    s = std(x(:,i));m = mean(x(:,i));
    x(:,i) = min(max(x(:,i),m-3*s),m+3*s);
    x(:,i) = x(:,i) - min(x(:,i));
    x(:,i) = ceil(nbins(min(i,numel(nbins)))*max(eps,x(:,i))/max(x(:,i)));
    x(:,i) = min(nbins(min(i,numel(nbins))),x(:,i));
end
data = accumarray(x,y,nbins,@mean);
al = accumarray(x,y,nbins,@ste);
al = accumarray(x,y,nbins,@histo);
al = al-min(al(:));
al = al/max(al(:));
data(data == 0) =nan;

function a = histo(a)
a = numel(a);

function a = gsorth(a)
a=complex(real(a),imag(a)-real(a)*diag(sum(real(a).*imag(a))./sum(real(a).^2)));

function x1 = ste(x)
x1 = (numel(x))/var(x);
x1 = log(x1);
if numel(x) < 2
    x1 = 0;
elseif abs(x1) == inf || isnan(x1)
    x1 = 0;
end