function [spT spId] = spikePlots(file,pos,v,spT,spId)%,sp)
%% make plots of spiking during linear track runs
fs = 1250/32;
bounds = [.1 .9];
accumbins = 50;timeBins = [-100:400];
angCol = colormap('hsv');
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:2
    pos(:,i) = interp1(find(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);%sp = sp(:,~nanInds);
nanInds = find(nanInds);
if ~exist('spT','var')
    d = '/media/Expansion Drive/KenjiMizuseki/';%['/media/work/hippocampus/' file '/'];%
    [spT,spId,~,d] = LoadCluRes([d file]);
    spT = (spT/d.SampleRate*fs);
    spT = spT - max(nanInds(nanInds < size(pos,1)/2));
    inds = spT < 1 | spT > max(size(pos,1)); spT(inds) = [];spId(inds) = [];
end

vel = angVel(pos);vel = filtLow(vel(:,1),fs,1);
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
figure;
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    temp = spId(ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0)));
    [~,spSort] = sort(hist(temp,1:max(spId)),'descend');
for j = 1:max(spId)
%for i = 1:max(runs)
    inds = find(runs == i);inds = min(inds):max(inds);
    %hold off;plot(0,0);
    %for j = 18%:20%1:max(spId)
 %       tempTimes = spT(spT > min(inds) & spT < max(inds) & spId == spSort(j));
 tempTimes = spT(spId == spSort(j));tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
        %scatter(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,v(:,1))));drawnow;pause(1)
        %scatter(weighted(tempTimes,pos(:,1)),weighted(tempTimes,(1:size(pos,1))'),30,angCol(phase2Col(angle(weighted(tempTimes,v(:,1)))),:),'filled');
        %plot(hist(weighted(tempTimes,pos(:,1)),linspace(bounds(1),bounds(2),30)));hold all;
        imagesc(log(hist3([angle(weighted(tempTimes,v(:,1))) weighted(tempTimes,pos(:,1))],[50 50])));
        title(num2str(numel(tempTimes)));colorbar;pause(.5);
        %temp = pos(sp(spSort(j),:) >0 & runs > 0,1);
        %hold all;plot(hist(temp,linspace(bounds(1),bounds(2),30)),'r');hold off;
%        hold all;
%    end
%    drawnow;pause(.1);hold off;scatter(0,0);
    %indsa = min(inds)-100:max(inds);
    %start = find(vel(indsa) > thresh,1);
    %start = max(min(indsa)+start-1,-min(timeBins)+1);
    %inds(vel(inds,1) < .1) = [];
end
end

function c = phase2Col(ang)
c = ceil((ang+pi)/(2*pi)*64);

function out = weighted(inds,vals)
w = mod(inds,1);
out = (1-w).*vals(floor(inds))+w.*vals(ceil(inds));