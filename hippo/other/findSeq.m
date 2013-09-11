function seq = findSeq(allCorr,ind,thresh)
%%an attempt to track sequence of activations in waiting area for non
%%place-locked FFP's.

allCorr = abs(allCorr);
mid = (size(allCorr,1)+1)/2;
multVec = 1;%[zeros(1,mid-1) ones(1,mid)];%1 exp(-(0:mid-2)/mid)];%
allCorr(:,:) = bsxfun(@times,allCorr(:,:),multVec');
meth = 1;
if meth == 0
allCorr = max(0,allCorr - thresh(1));
cumCorr = squeeze(sum(allCorr));
cumCorr = cumCorr;%cumCorr - thresh);
cumTime = zeros(size(cumCorr));
for i = 1:size(allCorr,1)
    cumTime = cumTime + (i-(size(allCorr,1)+1)/2)*squeeze(allCorr(i,:,:));
end
cumTime = cumTime./cumCorr;
else
[cumCorr,cumTime] = max(allCorr);
cumCorr = squeeze(cumCorr);
cumTime = squeeze(cumTime)-(size(allCorr,1)+1)/2;
cumCorr = max(0,cumCorr - thresh);
for i = 1:size(cumCorr,1)
    cumCorr(i,i) = 0;
    temp = sort(cumCorr(i,:),'descend');
    cumCorr(i,cumCorr(i,:) < temp(2)) = 0;
end
thresh = 0;
end
%cumCorr(cumTime < 0) = 0;
oldInds = [];
newInds = ind;
seqTimes = nan*ones(size(cumCorr,1),1);
seqTimes(ind) = 0;
counter = 0;
figure;scatter(0,0,'filled');hold all;drawnow;
while numel(newInds)
    counter = counter + 1;
    addInds = [];
    for i= 1:numel(newInds)
%%        [t,temp] = sort(cumCorr(:,newInds(i)),'descend');
%%        if t(1) > thresh(end)
%%            addInds = [addInds temp(1)];
%%        end
        addInds = [addInds find(cumCorr(:,newInds(i)) > thresh(end))'];
    end
    addInds = unique(addInds)
    oldInds = union(newInds,oldInds);
%    [fx,fy] = find(cumCorr(:,oldInds) == max(max(cumCorr(:,oldInds))));
%    if max(max(cumCorr(:,oldInds))) > thresh
%        addInds = fx;
%    end
    newInds = addInds(~ismember(addInds,oldInds));
    for i = 1:numel(newInds)
        refInds = intersect(oldInds,find(cumCorr(newInds(i),:) > thresh(end)));
        seqTimes(newInds(i)) = sum((seqTimes(refInds)'+cumTime(newInds(i),refInds))...
            .*cumCorr(newInds(i),refInds))/sum(cumCorr(newInds(i),refInds));
    end
%    if numel(newInds)
%        cumCorr(fx,:) = 0;
%    end
    scatter(counter*ones(1,numel(newInds)),seqTimes(newInds),'filled');drawnow;
end
[~,temp] = sort(seqTimes(oldInds),'ascend');
seq = oldInds(temp);%temp(~isnan(seqTimes));
%seq = nan*ones(size(seqTimes));