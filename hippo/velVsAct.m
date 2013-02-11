function [actVal velVal] = velVsAct(act,vel,cull,pad)

if cull > 0
    for i = 1:size(act,1)
        temp = squeeze(act(i,:,:));
        t = mean(temp);
        aa = bsxfun(@minus,temp,t);%temp - u*s*v1';
        snr(i) = sum(abs(temp(:)).^2)./sum(abs(aa(:)).^2);
    end
    snr(isnan(snr)) = 0;
    [~,inds] = sort(snr,'descend');
    act = act(inds(1:cull),:,:);
end
% figure;plot(snr);
% figure;
% for i = 1:size(act,1)
%     subplot(6,6,i);imagesc(squeeze(abs(act(i,:,:))));
% end
actm = abs(squeeze(mean(act,2)));
[~,maxInds] = max(actm,[],2);
figure;
for i = 1:size(act,1)
    %temp = squeeze(act(i,:,:));
    weighted = actm(i,:);
    subplot(211);plot(weighted);hold all;
    weighted(1:maxInds(i)-pad) = 0;
    weighted(maxInds(i)+pad:end) = 0;
    subplot(212);plot(weighted);hold all;
    actVal(i,:) = sum(bsxfun(@times,squeeze(abs(act(i,:,:))),weighted),2);
    actVal(i,:) = actVal(i,:)/sum(weighted);
    velVal(i,:) = sum(bsxfun(@times,vel,weighted),2);
    velVal(i,:) = velVal(i,:)/sum(weighted);
end
