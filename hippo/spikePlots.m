function [spT spId spf] = spikePlots(file,pos,v,spT,spId,spf)%,sp)
%% make plots of spiking during linear track runs
fs = 1250/32;
bounds = [.1 .9];
accumbins = 50;timeBins = [-100:400];
%angCol = colormap('hsv');
pos(pos == -1) = nan;
sV = size(v,1);
%if size(v,1) < size(pos,1)
pos = pos(1:sV,:);
%end
for i = 1:2
    pos(:,i) = interp1(find(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);
if ~exist('spT','var')
    %d = ['/media/Expansion Drive/KenjiMizuseki/'; file '/'];%%'/media/work/hippocampus/'
    %[spT,spId,~,d] = LoadCluRes([d file],[1:4 7:10]);
    %spT = (spT/d.SampleRate*fs);
    if strcmp('ec013.670',file)
        shanks = 5:8;
    elseif strcmp('ec014.468',file)
        shanks = 1:8;
    elseif strcmp('ec016.269',file)
        shanks = [1:4 7:10];
    end
    [spf spT spId] = hipSpikes(file,1000/fs,shanks);
    spf = morFilter(spf,8,fs);
    return;
end
spf = spf(:,1:sV);
spf = spf(:,~nanInds);
nanInds = find(nanInds);
spT = spT - max(nanInds(nanInds < size(pos,1)/2));
inds = spT < 1 | spT > max(size(pos,1)); spT(inds) = [];spId(inds) = [];
spf = bsxfun(@times,spf,exp(1i*angle(v(:,1))).');
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
vDemod(:,1) = [0; v(2:end,1).*exp(1i*-angle(v(1:end-1,1)))];
vDemod(:,2) = v(:,2).*exp(1i*-angle(v(:,1)));
bins{1} = linspace(-pi,pi,30);bins{2} = bins{1};%
binsp{1} = bins{1};binsp{2} = linspace(bounds(1),bounds(2),50);%
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    temp = spId(ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0)));
    [~,spSort] = sort(hist(temp,1:max(spId)),'descend');
    posHist = hist(pos(runs > 0,1),binsp{2});
    v2Hist = hist(angle(vDemod(runs > 0,2)),bins{2});
    %h0 = hist3([angle(vDemod(runs>0,2)) pos(runs>0,1)],binsp);
    %h0 = bsxfun(@rdivide,h0,posHist);
    %h1 = hist3(angle(vDemod(runs>0,:)),bins);
    %h1 = bsxfun(@rdivide,h1,v2Hist);
for j = 1:max(spId)
%for i = 1:max(runs)
%    inds = find(runs == i);inds = min(inds):max(inds);
 tempTimes = spT(spId == spSort(j));tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
 if numel(tempTimes) > 100
       hp = hist3([angle(weighted(tempTimes,spf(spSort(j),:).')) weighted(tempTimes,pos(:,1))],binsp);
       hp = bsxfun(@rdivide,hp,posHist);
       hp2 = hist3([angle(weighted(tempTimes,(spf(spSort(j),:).*exp(1i*angle(vDemod(:,2).'))).')) weighted(tempTimes,pos(:,1))],binsp);
       hp2 = bsxfun(@rdivide,hp2,posHist);
       hv2p = hist3([angle(weighted(tempTimes,v(:,2))) weighted(tempTimes,pos(:,1))],binsp);
       hv2p = bsxfun(@rdivide,hv2p,posHist);
       h = hist3([angle(weighted(tempTimes,v(:,1))) angle(weighted(tempTimes,vDemod(:,2)))],bins);
        figure(2);subplot(221);imagesc(sqrt(hp));subplot(222);imagesc(sqrt(hp2));
        subplot(223);imagesc(sqrt(h));
        h = bsxfun(@rdivide,h,v2Hist);
        subplot(224);imagesc(sqrt(h));
%        subplot(224);imagesc(sqrt(hv2p));
         input('7');
drawnow;
 end
end
axis tight;
end

function c = phase2Col(ang)
c = ceil((ang+pi)/(2*pi)*64);

function out = weighted(inds,vals)
w = mod(inds,1);
out = (1-w).*vals(floor(inds))+w.*vals(ceil(inds));