function [spT spId spf] = spikePlots(file,pos,v,spT,spId,spf)%,sp)
%% make plots of spiking during linear track runs
fs = 1250/32;
bounds = [.1 .9];
accumbins = 50;timeBins = [-100:400];
angCol = colormap('hsv');
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
    [spf spT spId] = hipSpikes(file,1000/fs,[5 6 11 12]);
    spf = morFilter(spf,8,fs);
    size(spf)
    sV
    spf = spf(:,1:sV);
    spf = spf(:,~nanInds);
    nanInds = find(nanInds);
    spT = spT - max(nanInds(nanInds < size(pos,1)/2));
    inds = spT < 1 | spT > max(size(pos,1)); spT(inds) = [];spId(inds) = [];
    spf = bsxfun(@times,spf,exp(1i*angle(v(:,1))).');
end

%vel = angVel(pos);vel = filtLow(vel(:,1),fs,1);
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
v(:,2) = v(:,2).*exp(1i*-angle(v(:,1)));
bins{1} = linspace(-pi,pi,30);bins{2} = linspace(bounds(1),bounds(2),50);%bins{1};%
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    temp = spId(ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0)));
    [~,spSort] = sort(hist(temp,1:max(spId)),'descend');
    posHist = hist(pos(runs > 0,1),bins{2});
    v2Hist = hist(angle(v(runs > 0,2)),bins{2});
    h0 = hist3([angle(v(runs>0,2)) pos(runs>0,1)],bins);
    h0 = bsxfun(@rdivide,h0,posHist);
for j = 1:max(spId)
%for i = 1:max(runs)
    inds = find(runs == i);inds = min(inds):max(inds);
    %hold off;plot(0,0);
    %for j = 18%:20%1:max(spId)
 %       tempTimes = spT(spT > min(inds) & spT < max(inds) & spId == spSort(j));
 tempTimes = spT(spId == spSort(j));tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
 if numel(tempTimes) > 100
        %scatter(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,v(:,1))));drawnow;pause(1)
%         figure(1);
%         %scatter(weighted(tempTimes,pos(:,1)),weighted(tempTimes,(1:size(pos,1))'),30,angCol(phase2Col(angle(weighted(tempTimes,spf(spSort(j),:).'))),:),'filled');
%         plot(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,v(:,1))),'b','linewidth',2);hold all;
%         plot(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,spf(spSort(j),:).')),'r','linewidth',2);
%         scatter(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,v(:,1))),'b','filled');hold all;
%         scatter(weighted(tempTimes,pos(:,1)),angle(weighted(tempTimes,spf(spSort(j),:).')),'r','filled');hold off;
%         %plot(hist(weighted(tempTimes,pos(:,1)),linspace(bounds(1),bounds(2),30)));hold all;
%         %plot(pos(inds,1),i/100+real(v(inds,1)),'b','linewidth',2);hold all;
%         %plot(pos(inds,1),i/100+imag(v(inds,1)),'r','linewidth',2);
       h = hist3([angle(weighted(tempTimes,spf(spSort(j),:).')) weighted(tempTimes,pos(:,1))],bins);
       h = bsxfun(@rdivide,h,posHist);
       h1 = hist3([angle(weighted(tempTimes,v(:,1))) weighted(tempTimes,pos(:,1))],bins);
       h1 = bsxfun(@rdivide,h1,posHist);
       %h = hist3([angle(weighted(tempTimes,v(:,1))) angle(weighted(tempTimes,v(:,2)))],bins);
       %h = bsxfun(@rdivide,h,v2Hist);
        figure(2);subplot(311);imagesc(sqrt(h0));
        subplot(312);imagesc(sqrt(h));
        subplot(313);imagesc(sqrt(h1));
%         %title(num2str(numel(tempTimes)));colorbar;pause(.5);
%         %temp = pos(sp(spSort(j),:) >0 & runs > 0,1);
%         %hold all;plot(hist(temp,linspace(bounds(1),bounds(2),30)),'r');hold off;
         input('7');%
%        pause(1);
drawnow;
    end
%    drawnow;pause(.1);hold off;scatter(0,0);
    %indsa = min(inds)-100:max(inds);
    %start = find(vel(indsa) > thresh,1);
    %start = max(min(indsa)+start-1,-min(timeBins)+1);
    %inds(vel(inds,1) < .1) = [];
% end
end
axis tight;
end

function c = phase2Col(ang)
c = ceil((ang+pi)/(2*pi)*64);

function out = weighted(inds,vals)
w = mod(inds,1);
out = (1-w).*vals(floor(inds))+w.*vals(ceil(inds));