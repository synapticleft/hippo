function [spf spT spId spShank] = spikePlots(file,pos,v,spf,spT,spId)%,sp)
%% make plots of spiking during linear track runs
fs = 1250/32;               %Sampling Rate
bounds = [.1 .9];           %Part of track to consider runway
accumbins = 50;             %binning of track length
pos(pos == -1) = nan;
%% Discretize position along track
sV = size(v,1);
%if size(v,1) < size(pos,1)
pos = pos(1:sV,:);
%end
for i = 1:2
    pos(:,i) = interp1(find(~isnan(pos(:,i))),pos(~isnan(pos(:,i)),i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/(max(pos(:,1))+.00001);
posd = floor(pos(:,1)*accumbins)+1;
%% extract spikes from file if havent done so already
if ~exist('spf','var')
    if numel(strfind(file,'ec013'))%strcmp('ec013.670',file)
        shanks = 5:8;
    elseif numel(strfind(file,'ec014'))%strcmp('ec014.468',file)
        shanks = 1:8;
    elseif numel(strfind(file,'ec016'))%strcmp('ec016.269',file)
        shanks = [1:4 7:10];
    end
    [spf spT spId spShank bad] = hipSpikes(file,1000/fs,shanks);
    spf = [bad;spf];            %lump artifacts & unsorted spikes at beginning of matrix
    spf = morFilter(spf,8,fs);  %hilbert-transform spiking data
    return;
end
spf = spf(:,1:sV);
spf = spf(:,~nanInds);
nanInds = find(nanInds);
if exist('spT','var')
    spT = spT - max(nanInds(nanInds < size(pos,1)/2));
    inds = spT < 1 | spT > max(size(pos,1)); spT(inds) = [];spId(inds) = [];
end
%% Demodulate filtered spike train
spf = bsxfun(@times,spf,exp(1i*angle(v(:,1))).');
%% separate runs into trials
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
runs = bwlabel(b > 0);
spInterp = zeros(2,size(spf,1),max(runs),accumbins(1));
w = watershed(b==0);
w = w-1; %w(w== max(w)) = 0;
%% accummulate spiking samples for every trial and position bin
for k = 1:2
    runs = bwlabel(w>0 & mod(w,2) == k-1);%b*((-1)^k)>0);
    inds = runs > 0;
    for j = 1:size(spf,1)
        spInterp(k,j,:,:) = accumarray([runs(inds); posd(inds)']',spf(j,inds),[max(runs) accumbins] ,@mean);
    end
end
spInterp = [squeeze(spInterp(1,:,:)) squeeze(spInterp(2,:,:))];
%% display results
ydim = 8;xdim = ceil(size(spf,1)/ydim);figure;
for i = 1:size(spInterp,1)
    temp = reshape(spInterp(i,:),[max(runs) 2*accumbins]);
    subplot(xdim,ydim,i);imagesc(complexIm(temp,0,1));axis off;
end

%% introns
% vDemod(:,1) = [0; v(2:end,1).*exp(1i*-angle(v(1:end-1,1)))];
% vDemod(:,2) = v(:,2).*exp(1i*-angle(v(:,1)));
% angBins = 20;posBins = 30;
% bins{1} = linspace(-pi+pi/angBins,pi-pi/angBins,angBins);bins{2} = bins{1};%
% binsp{1} = bins{1};binsp{2} = linspace(bounds(1)+1/posBins/2,bounds(2)-1/posBins/2,posBins);%
% binspab = binsp;binspab{1} = linspace(0,2*std(vDemod(:,2)),angBins);
% for k = 1:2
%     runs = bwlabel(b*((-1)^k)>0);
%     posHist = hist(pos(runs > 0,1),binsp{2});
%     mem = ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0));
%     subTimes = spT(mem);subId = spId(mem);
%     [~,spSort] = sort(sum(abs(spf(:,runs > 0)),2),'descend');
% %     runs = bwlabel(b*((-1)^k)>0);
% %     temp = spId(ismember(floor(spT),find(runs>0)) | ismember(ceil(spT),find(runs>0)));
% %     [~,spSort] = sort(hist(temp,1:max(spId)),'descend');
% %     posHist = hist(pos(runs > 0,1),binsp{2});
%     v2Hist = hist(angle(vDemod(runs > 0,2)),bins{2});
%     h0 = hist3([angle(vDemod(runs>0,2)) pos(runs>0,1)],binsp);
%     h0 = bsxfun(@rdivide,h0,posHist);
%     h1 = hist3([angle(v(runs>0,1)) angle(vDemod(runs>0,2))],bins);%angle(vDemod(runs>0,:)),bins);
%     %h1 = bsxfun(@rdivide,h1,v2Hist);
%     h2 = hist3([abs(vDemod(runs>0,2)) pos(runs>0,1)],binspab);
%     h2 = bsxfun(@rdivide,h2,posHist);
%     figure(1);subplot(2,2,1);imagesc(h0);subplot(2,2,2);imagesc(h1);
%     subplot(2,2,3);imagesc(h2);
%     h3 = hist3([angle(v(runs>0,:))],bins);
%     subplot(2,2,4);imagesc(h3);
% for j = 1:10%max(spId)
% %for i = 1:max(runs)
% %    inds = find(runs == i);inds = min(inds):max(inds);
%  tempTimes = spT(spId == spSort(j));tempTimes = tempTimes(ismember(floor(tempTimes),find(runs>0)) | ismember(ceil(tempTimes),find(runs>0)));
%  if numel(tempTimes) > 100
%        hp = hist3([angle(weighted(tempTimes,spf(spSort(j),:).')) weighted(tempTimes,pos(:,1))],binsp);
%        hp = bsxfun(@rdivide,hp,posHist);
%        hp2 = hist3([angle(weighted(tempTimes,(spf(spSort(j),:).*exp(1i*angle(vDemod(:,2).'))).')) weighted(tempTimes,pos(:,1))],binsp);
%        hp2 = bsxfun(@rdivide,hp2,posHist);
%        hv2p = hist3([angle(weighted(tempTimes,v(:,2))) weighted(tempTimes,pos(:,1))],binsp);
%        hv2p = bsxfun(@rdivide,hv2p,posHist);
%        h = hist3([angle(weighted(tempTimes,v(:,1))) angle(weighted(tempTimes,vDemod(:,2)))],bins);
%         figure(2);subplot(231);imagesc((hp));subplot(232);imagesc((hp2));
%         temp = zeros(size(hp));
%         for i = 1:size(temp,2)
%             temp(:,i) = cconv(hp(:,i),h0(:,i),size(temp,1));
%         end
%         subplot(233);imagesc(circshift(temp,[size(temp,1)/2 0]));
%         subplot(234);imagesc((h));
%         subplot(236);plot(mean(h)/mean(h(:)));hold all;plot(v2Hist/mean(v2Hist));hold off;%scatter(mean(h),v2Hist);
%         h = bsxfun(@rdivide,h,v2Hist);
%         subplot(235);imagesc((h));
% %        subplot(224);imagesc(sqrt(hv2p));
%          input(num2str(spSort(j)));
% drawnow;
%  end
% end
% axis tight;
% end

% function out = weighted(inds,vals)
% w = mod(inds,1);
% out = (1-w).*vals(floor(inds))+w.*vals(ceil(inds));