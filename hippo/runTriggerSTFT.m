function runTriggerSTFT(pos,X,thresh)
%downsample position vector so that can look at phase relations of theta in
%v
warning off all;
bounds = [.1 .9];
accumbins = 50;timeBins = -10:40;
pos(pos == -1) = nan;
% if size(v,1) < size(pos,1)
%     pos = pos(1:size(v,1),:);
% end
nanInds = find(~isnan(pos(:,1)));
for i = 1:4
    nans = find(~isnan(pos(:,i)));
pos(:,i) = interp1(nans,pos(nans,i),1:size(pos,1));
end
%pos(:,2) = interp1(nanInds,pos(nanInds,2),1:size(pos,1));
nanInds = isnan(pos(:,1));
%%STFT
params.Fs = 1250/25;
if exist('dec','var')
    params.Fs = params.Fs/dec;
end
win = [.5 .3];
params.tapers = [1/win(1) win(1) 1];
[S,t,f,sa] = mtspecgramc(X,win,params);
ref = 5:7;ref1 = 12:13;
ft = exp(-1i*angle(mean(sa(:,ref),2))*ones(size(f)));
temp = sa.*ft;
%%resample
pos(nanInds,:) = 0;
pos = resample(pos,size(temp,1),size(pos,1));
nanInds = resample(double(nanInds),size(temp,1),size(nanInds,1));
pos = pos(nanInds < 0.1,:);
temp = temp(nanInds < 0.1,:);
sa = sa(nanInds < .1,:);
vel = angVel(pos);vel = filtLow(vel(:,1),1/win(2),1);
figure;image(complexIm(temp.',0,1));hold all;plot(vel/10,'linewidth',2);
pos = bsxfun(@minus,pos,mean(pos));%pos = bsxfun(@rdivide,pos,std(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
v = [mean(sa(:,ref),2) mean(sa(:,ref(2)),2)];%12:13
v(:,2) = v(:,2).*(conj(v(:,1))./abs(v(:,1)));%.^2;
v(:,1) = [0; v(2:end,1).*conj(v(1:end-1,1))./abs(v(1:end-1,1))];
if 1/(2*win(2)) > 1
v = filtLow(v.',1/win(2),1).';
end
runs = bwlabel(b > 0);
vInterp = zeros(2,2,max(runs),accumbins);
velInterp = zeros(2,max(runs),accumbins);
velTrace = zeros(2,max(runs),range(timeBins)+1);
vVel = zeros(2,2,max(runs),range(timeBins)+1);
bins = (bounds(1))+((1:accumbins)-.5)/accumbins*(diff(bounds));
%plot(pos(:,1));hold all;plot(vel);return
%figure;plot(vel/10);hold all;plot(abs(b));
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    %avg = zeros(size(sp,1),accumbins,max(lrRuns));
for i = 1:max(runs)
    inds = find(runs == i);inds = min(inds):max(inds);
    indsa = min(inds)+timeBins(1):max(inds);
    start = find(vel(indsa) > thresh,1);
    if isempty(start)
        start = -timeBins(1);
    end
    start = max(min(indsa)+start-1,-min(timeBins)+1);
    velTrace(k,i,:) = vel(start+timeBins);
    inds(vel(inds,1) < .1) = [];
    %for j = 1:size(sp,1)
    %    lrAvg(j,:,i) = csaps(pos(inds,1),sp(j,inds),1-1e-7,bins);
    %end
    for j = 1:2
        vVel(k,j,i,:) = v(start+timeBins,j);
        vInterp(k,j,i,:) = csaps(pos(inds,1),v(inds,j),1-1e-7,bins);
    end
    velInterp(k,i,:) = csaps(pos(inds,1),vel(inds),1-1e-7,bins);
end
end
numGraphs = 6;
figure;
powers = [1 1];
%vInterp(:,:,2,1:6) =0;
for i = 1:2
    subplot(2,numGraphs,(i-1)*numGraphs+1);imagesc(timeBins*32/1250,1:max(runs),squeeze(velTrace(i,:,:)),[0 prctile(velTrace(:),99)]);
    subplot(2,numGraphs,(i-1)*numGraphs+6);imagesc(linspace(bounds(1),bounds(2),accumbins)*250,1:max(runs),squeeze(velInterp(i,:,:)),[0 prctile(velInterp(:),99)]);
    set(gca,'fontsize',16);title 'velocity';ylabel 'Trial #'; xlabel 'Time (ms)'
    for j = 1:2
        subplot(2,numGraphs,(i-1)*numGraphs+2*j);imagesc(timeBins*32/1250,1:max(runs),complexIm(squeeze(vVel(i,3-j,:,:)),0,powers(j)));
        set(gca,'fontsize',16);
        subplot(2,numGraphs,(i-1)*numGraphs+2*j+1);imagesc(linspace(bounds(1),bounds(2),accumbins)*250,1:max(runs),complexIm(squeeze(vInterp(i,3-j,:,:)),0,powers(j)));
        set(gca,'fontsize',16);
    end
end