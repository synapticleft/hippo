function posVsTime(pos,v,bounds,inds)%,sp

accumbins = [60 10];
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
vp12 = vp11;%v(:,1).*conj(v(:,2))./abs(v(:,1));%
state = [pos(b>0); 1:numel(pos(b>0))]';
[datalr ] = binData(state,vp12(b>0),accumbins);
state = [pos(b<0); 1:numel(pos(b<0))]';
[datarl ] = binData(state,vp12(b<0),accumbins);
figure;subplot(211);imagesc(abs(datalr)');
subplot(212);imagesc(abs(datarl)');
figure;subplot(211);imagesc(angle(datalr)');
subplot(212);imagesc(angle(datarl)');colormap hsv
%%%%%
statelr = [pos(b>0)' angle(vp12(b>0))];
staterl = [pos(b<0)' angle(vp12(b<0))];
numBins= [30 50];numDivs = 1;
makeStateFig(statelr,staterl,numBins,numDivs);
%statelr = [pos(b>0)' abs(vp12(b>0))];
%staterl = [pos(b<0)' abs(vp12(b<0))];
%makeStateFig(statelr,staterl,numBins,numDivs);

function makeStateFig(statelr,staterl,numBins,numDivs)
figure;
for i = 1:numDivs
    inds = floor((i-1)*size(statelr,1)/numDivs+1):ceil(i*size(statelr,1)/numDivs);
    [h bi] = hist3(statelr(inds,:),numBins);
	subplot(2,numDivs,i);imagesc(bi{1},bi{2},h');
    inds = floor((i-1)*size(staterl,1)/numDivs+1):ceil(i*size(staterl,1)/numDivs);
    [h bi] = hist3(staterl(inds,:),numBins);
    subplot(2,numDivs,numDivs+i);imagesc(bi{1},bi{2},h');
end
drawnow;

function [data al] = binData(x,y,nbins)
for i = 1:size(x,2)
    s = std(x(:,i));m = mean(x(:,i));
    x(:,i) = min(max(x(:,i),m-2*s),m+2*s);
    x(:,i) = x(:,i) - min(x(:,i));
    x(:,i) = ceil(nbins(min(i,numel(nbins)))*max(eps,x(:,i))/max(x(:,i)));
    x(:,i) = min(nbins(min(i,numel(nbins))),x(:,i));
end
data = accumarray(x,y,nbins,@mean);
al = accumarray(x,y,nbins,@std);
data(data == 0) =nan;