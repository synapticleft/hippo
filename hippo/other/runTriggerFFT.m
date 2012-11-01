function runTriggerFFT(FT,f,pos)

[~,ind] = min(abs(f-8));
[u,~,v] = svds(squeeze(FT(:,ind,:)).',1);
FT = realign(FT,f,v,ind);
figure;imagesc(f,f,complexIm(squeeze(mean(FT,3)),0,1));

function ft = realign(ft,f,v,ind)
for i = 1:numel(f)
    scale = f(i)/f(ind);
    ft(:,i,:) = bsxfun(@times,squeeze(ft(:,i,:)),exp(1i*angle(v)*scale));
end