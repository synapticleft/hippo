function restlessClassify(act,pos,bounds)
%% Comparing waiting area decoders for discerning previous and future choices
%% Ans: previous easier to decode

thresh = 3;
dec = round(size(pos,1)/size(act,2));
for i =1:2
    posd(:,i) = decimate(pos(:,i),dec);
end
posd = posd(1:size(act,2),:);
bounded = posd(:,1) > bounds(1) & posd(:,1) < bounds(2) & ...
    posd(:,2) > bounds(3) & posd(:,2) < bounds(4);
numExp = 13*8/dec;
bounded = bwmorph(bounded,'dilate',numExp);
for i =1:numExp
    bounded = bounded - bwmorph(bounded,'endpoints');
end
bounded = logical(bounded);
bounded = bwlabel(bounded);
h = hist(bounded,0:max(bounded));h(1) = [];bounded(ismember(bounded,find(h < 10))) = 0;
bounded = bwlabel(bounded > 0);
[posMean posSeq times] = sideSeq(bounded,posd,mean(bounds(3:4)));
%posSeq = posSeq(randperm(numel(posSeq)));
%subplot(122);scatter(posd(bounded > 0,1),posd(bounded > 0,2),'b','filled');hold on;
figure;hold on;
val = zeros(size(act,1),max(bounded));
fs = find(times(1:end-1) > 30);
err= zeros(size(act,1),1);
err1 = err;
for i = 1:size(act,1)
    val(i,:) = (accumarray(bounded(bounded > 0),abs(act(i,bounded > 0)),[max(bounded) 1],@max,nan));
    [~,err(i)] = classify(val(i,:)',val(i,fs)',posSeq(fs)>0);
    [~,err1(i)] = classify(val(i,:)',val(i,fs)',posSeq(fs+1)>0);
    %valShift = circshift(val,[1 0]);
    %subplot(131);scatter(log(times),val,'filled');hold all;
    %scatter(log(times(posSeq > 0)),val(posSeq > 0),'r','filled');hold off;
    %subplot(132);scatter(log(times),valShift,'filled');hold all;
    %scatter(log(times(posSeq > 0)),valShift(posSeq > 0),'r','filled');hold off;
    %subplot(133);scatter(posd(bounded > 0,1),posd(bounded > 0,2),'b','filled');hold all;
    %s = scatter(posd(bounded > 0 & abs(act(i,:)') > thresh,1),...
    %    posd(bounded > 0 & abs(act(i,:)') > thresh,2),'r','filled');hold off;drawnow;input('');
end
scatter(err,err1,'b','filled');
[c1,err] = classify(val',val(:,fs)',posSeq(fs)>0,'diaglinear');
[c2,err1] = classify(val',val(:,fs)',posSeq(fs+1)>0,'diaglinear');
scatter(err,err1,'r','filled');hold on;plot([0 .7],[0 .7],'k--');
set(gca,'fontsize',16);title '% classifier error';
xlabel('On previous choice');
ylabel('On future choice');
legend({'Individual TMPCs';'All TMPCs'});axis square tight;
figure;plot(posd(:,1:2));hold all;
%plot((bounded > 0)*300);
pos1 = nan*ones(size(bounded));
pos2 = pos1;
for i = 1:max(bounded)
    pos1(bounded == i) = c1(i);%posSeq(i)>0;%
    pos2(bounded == i) = c2(i);%posSeq(i+1)>0;%
end
plot(pos1*50+200+randn(size(pos1)));
plot(pos2*50+200+randn(size(pos1)));
%plot((bounded > 0)*300);
% for i = 1:2
%     figure;hold all;
%     inds = bounded & posMean == i-1;
%     for j = 1:size(act,1)
%     [ys xs] = meshgrid(find(inds),j);
%     tempMag = abs(act(j,inds));
%     scatter(posd(ys(tempMag > thresh),2),posd(ys(tempMag > thresh),1),...
%         tempMag(tempMag > thresh)*3,getCol(xs,tempMag>thresh,frameCol),'filled');
%     drawnow;
%     end
%     set(gca,'xlim',max(0,bounds(3:4)),'ylim',bounds(1:2),'color',[0 0 0]);
%     axis square
% end

function [runDir posSeq times] = sideSeq(runInds,pos,ref)
posSeq = zeros(1,max(runInds)+1);
runDir = zeros(size(runInds));
times = posSeq;
for i = 1:max(runInds)+1
    if i == 1
        inds = 1:find(runInds == 1,1)-1;
    elseif i == max(runInds)+1
        inds = find(runInds == i-1, 1, 'last' )+1:size(pos,1);
    else
        inds = (find(runInds == i-1,1,'last')+1):(find(runInds == i,1)-1);
    end
    temp = [max(pos(inds,2))-ref min(pos(inds,2))-ref];
    times(i) = sum(runInds == i);
    if abs(temp(1)) > abs(temp(2))
        posSeq(i) = temp(1);
    else
        posSeq(i) = temp(2);
    end
    runDir(runInds == i) = posSeq(i) > 0;
end