function stateHists(pos,vars)

if size(pos,1) < numel(vars)
    vars((size(pos,1)+1):end) = [];
else
    pos((numel(vars)+1):end,:) = [];
end
t = (1:size(pos,1))*32/1250;
win = 200;
pos(pos == -1) = nan;
for i = 1:size(pos,2)
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = any(isnan(pos'));
pos = pos(~nanInds,:);
pos = filtfilt(gausswin(win),1,pos);
pos = bsxfun(@minus,pos,mean(pos));
pos = bsxfun(@rdivide,pos,std(pos));
[u s v] = svd(pos(:,[1 3]),'econ');
%figure;plot(binData(pos(:,1),[diff(pos(:,1)); 0],30,u(:,2) > 0));
%figure;imagesc(hist3(pos(:,[1 2]),[10 30]))
if exist('vars','var')
    vars = vars(~nanInds);
    figure;plot(binData(pos(:,1),vars,30,u(:,2) > 0));
end

function data = binData(x,y,nbins,path)
x = x - min(x);
x = ceil(nbins*max(eps,x)/max(x));
data(:,1) = accumarray(x,y,[],@mean);
if exist('path','var')
    data(:,2) = accumarray(x(path),y(path),[],@mean);
    data(:,3) = accumarray(x(~path),y(~path),[],@mean);
end