function [absMean] = runTriggerView2dHigh(pos,Xf,scale,binSize)
%% binning of 2-d data in open maze

pos = fixPos(pos);
pos = pos(:,1:2);
for i = 1:2
    pos1(:,i) = resample(pos(:,i),scale,1);
end
pos = pos1(1:size(Xf,2),:);
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i))+eps;
    posd(:,i) = ceil(pos(:,i)/binSize);
end
absMean = zeros(size(Xf,1),max(posd(:,1)),max(posd(:,2)));
for i = 1:size(Xf,1)
    absMean(i,:,:) = accumarray(posd,(Xf(i,:)),[],@mean);
end