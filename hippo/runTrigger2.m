function runTrigger2(pos,v,thresh)
bounds = [.1 .9];
accumbins = 50;timeBins = [-100:400];
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
nanInds = find(~isnan(pos(:,1)));
pos(:,1) = interp1(nanInds,pos(nanInds,1),1:size(pos,1));
pos(:,2) = interp1(nanInds,pos(nanInds,2),1:size(pos,1));
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);%sp = sp(:,~nanInds);
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);
pos = bsxfun(@minus,pos,mean(pos));%pos = bsxfun(@rdivide,pos,std(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
v(:,2) = v(:,2).*conj(v(:,1))./abs(v(:,1));
v(:,1) = [0; v(2:end,1).*conj(v(1:end-1,1))./abs(v(1:end-1,1))];
v = filtLow(v.',1250/32,4).';
runs = bwlabel(b > 0);
vInterp = zeros(2,2,max(runs),accumbins);
velTrace = zeros(2,max(runs),range(timeBins)+1);
vVel = zeros(2,2,max(runs),range(timeBins)+1);
bins = (bounds(1))+((1:accumbins)-.5)/accumbins*(diff(bounds));
%plot(pos(:,1));hold all;plot(vel);return
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    %avg = zeros(size(sp,1),accumbins,max(lrRuns));
for i = 1:max(runs)
    inds = find(runs == i);inds = min(inds):max(inds);
    indsa = min(inds)-100:max(inds);
    start = find(vel(indsa) > thresh,1);
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
end
end
numGraphs = 5;
figure;
powers = [1 .25];
%vInterp(:,:,2,1:6) =0;
for i = 1:2
    subplot(2,numGraphs,(i-1)*numGraphs+1);imagesc(timeBins*32/1250,1:max(runs),squeeze(velTrace(i,:,:)),[0 prctile(velTrace(:),99)]);
    set(gca,'fontsize',16);title 'velocity';ylabel 'Trial #'; xlabel 'Time (ms)'
    for j = 1:2
        subplot(2,numGraphs,(i-1)*numGraphs+2*j);imagesc(timeBins*32/1250,1:max(runs),complexIm(squeeze(vVel(i,j,:,:)),0,powers(j)));
        set(gca,'fontsize',16);
        subplot(2,numGraphs,(i-1)*numGraphs+2*j+1);imagesc(linspace(bounds(1),bounds(2),accumbins)*250,1:max(runs),complexIm(squeeze(vInterp(i,j,:,:)),0,powers(j)));
        set(gca,'fontsize',16);
    end
end

   
%     
% vel = angVel(pos);
% vel = filtLow(vel(:,1),1250/32,2);
% vel(isnan(vel)) = 0;
% plot(vel);hold all;
% vel = toeplitz(vel,nan*ones(back,1));
% x = find(vel(:,1) > thresh & ~any(vel(:,2:end)' > thresh)');
% %plot(vel(:,1));hold all;scatter(x,vel(x,1),'r');
% vel = vel(:,1);
% vAll = zeros(numel(x),2*back+1,2);
% velAll = zeros(numel(x),2*back+1);
% for i = 1:numel(x)
%     vAll(i,:,:) = v(x(i)+(-back:back),:);
%     velAll(i,:) = vel(x(i)+(-back:back));
% end
% figure;imagesc(velAll);
% v1 = (vAll(:,:,1).*conj(vAll(:,:,1))./abs(vAll(:,:,1)));
% v11 = (vAll(:,2:end,1).*conj(vAll(:,1:end-1,1))./abs(vAll(:,1:end-1,1)));
% v12 = (vAll(:,:,2).*conj(vAll(:,:,1))./abs(vAll(:,:,1)));
% figure;subplot(221);image(complexIm(v11,0));
% subplot(222);image(complexIm(v11,1));
% subplot(223);image(complexIm(v12,0));
% subplot(224);image(complexIm(v12,1));
% figure;plot(mean(velAll));
% figure;plot(mean(v1));hold all;plot(imag(mean(v11)));plot(real(mean(v11)));
% plot(imag(mean(v12)));plot(real(mean(v12)));
%
%figure;imagesc(vel(x+floor(back/2),:));
%figure;plot(mean(vel(x+floor(back/2),:)));

