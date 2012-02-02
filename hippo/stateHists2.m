function stateHists2(pos,resp)
%2D histogram of neural response properties
%   Detailed explanation goes here
pos(pos == -1) = nan;
v = diff(pos);v = filtLow(v',1250/32,3)';
pos(1,:) = [];
state = [pos(:,1) v(:,1)];
%figure;imagesc(log(hist3(state,[100 100])));
[data al] = binData(state,resp,[40 40]);
h = imagesc(data);
%set(gca,'color','k');
%set(h,'AlphaData',al);

function [data al] = binData(x,y,nbins)
for i = 1:size(x,2)
    x(:,i) = x(:,i) - min(x(:,i));
    x(:,i) = ceil(nbins(min(i,numel(nbins)))*max(eps,x(:,i))/max(x(:,i)));
    x(:,i) = min(nbins(min(i,numel(nbins))),x(:,i));
end
data = accumarray(x,y,nbins,@mean);
al = accumarray(x,y,nbins,@ste);
al = al-min(al(:));
al = al/max(al(:));
data(data == 0) =nan;

function x1 = ste(x)
x1 = (numel(x))/var(x);
x1 = log(x1);
if numel(x) < 2
    x1 = 0;
elseif (abs(x) == inf) || isnan(x)
    x1 = 0;
end