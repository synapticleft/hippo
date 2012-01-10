function [pos vars] = stateHists(pos,vars)

if exist('vars','var')
    if size(pos,1) < numel(vars)
        vars((size(pos,1)+1):end) = [];
    else
        pos((numel(vars)+1):end,:) = [];
    end
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
params.Fs = 1250/32;%params.tapers = [8 15];
[S,f] = mtspectrumc(diff(pos(:,1)),params);figure;plot(f,log(S));
pos = filtLow(pos',params.Fs,4,8)';
%pos = filtfilt(gausswin(win),1,pos);
[S,f] = mtspectrumc(diff(pos(:,1)),params);hold all;plot(f,log(S));
pos = bsxfun(@minus,pos,mean(pos));
pos = bsxfun(@rdivide,pos,std(pos));
[u s v] = svd(pos(:,[1 3]),'econ');
%trial = cumsum(abs(diff(u(:,2) > 0)));
%trial(end)
figure;plot(binData(pos(:,1),[diff(pos(:,1)); 0],10,u(:,2) > 0));
figure;plot(abs(diff(u))*100);
%[u s v] = svd(pos(:,1:2),'econ');
%pos = u;%diff(u)*10;

%figure;imagesc(hist3(pos(:,[1 2]),[10 30]))
if exist('vars','var')
    vars = vars(~nanInds);
    hold all;plot(vars);
    figure;plot(binData(pos(:,1),vars,10,u(:,2) > 0));
end
figure;plot(pos(:,1)/std(pos(:,1)));hold all;plot(u(:,2)>0);
figure;plot(xcov(abs(diff([pos(:,1); 0])),abs(vars),1000,'coeff'));
hold all;plot(xcov(abs(diff(pos(:,1))),1000,'coeff'),'r');

function data = binData(x,y,nbins,path)
x = x - min(x);
x = ceil(nbins*max(eps,x)/max(x));
data(:,1) = accumarray(x,y,[],@mean);
if exist('path','var')
    data(:,2) = accumarray(x(path),y(path),[],@mean);
    data(:,3) = accumarray(x(~path),y(~path),[],@mean);
end