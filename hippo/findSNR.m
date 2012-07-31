function findSNR(sp,ic)

params.Fs = 1;
params.tapers = [3 5];

icNew = zeros(2,size(ic,1),size(ic,2),size(ic,3)/2);
icNew(1,:,:,:) = ic(:,:,1:size(ic,3)/2);
icNew(2,:,:,:) = ic(:,:,size(ic,3)/2+1:end);
n = numel(mtspectrumc(squeeze(sp(1,1,1,:)),params));
spSpec = zeros(2,2,size(sp,2),n);
icSpec = zeros(2,2,size(icNew,2),n);
for k = 1:2
    for i = 1:size(sp,2)
        sig = mean(squeeze(sp(k,i,:,:)));
        spSpec(1,k,i,:) = mtspectrumc(sig,params);
        spSpec(2,k,i,:) = mean(mtspectrumc(bsxfun(@minus,squeeze(sp(k,i,:,:)),sig)',params),2);
    end
    for i = 1:size(icNew,2)
        sig = mean(squeeze(icNew(k,i,:,:)));
        icSpec(1,k,i,:) = mtspectrumc(sig,params);
        icSpec(2,k,i,:) = mean(mtspectrumc(bsxfun(@minus,squeeze(icNew(k,i,:,:)),sig)',params),2);
    end
end

spSpec = squeeze(spSpec(:,1,:,:) + spSpec(:,2,:,:));
icSpec = squeeze(icSpec(:,1,:,:) + icSpec(:,2,:,:));

%figure;plot(squeeze(spSpec(1,:,:))');
%figure;plot(squeeze(icSpec(1,:,:))');
figure;plot(squeeze(spSpec(1,:,:)./spSpec(2,:,:))');
figure;plot(squeeze(icSpec(1,:,:)./icSpec(2,:,:))');