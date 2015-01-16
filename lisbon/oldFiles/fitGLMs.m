function fitGLMs(fn,lambda,prc)
% fit LDA and regularized linear regression of start- and end- aligned
% trials, look at performance, and weights.

file = dir('*.mat');
load(file(fn).name,'data');

info = reshape([data{2:end,[7 4 5]}],size(data,1)-1,[]);
info(:,3) = info(:,3) < 0;
class = 1;
inds = 6;

for i = 1:size(data,1)-1
    tot(i) = sum(~isnan(data{i+1,inds}));%
end
%tot = min(tot);
totpr = prctile(tot,prc);
fi = find(tot > totpr);
startAligned = zeros(numel(fi),totpr);
endAligned = startAligned;
for j = 1:numel(fi)%size(data,1)-1
    f = find(~isnan(data{fi(j)+1,inds}));
    startAligned(j,:) = data{fi(j)+1,inds}(f(1)+(0:totpr-1));
    endAligned(j,:) = data{fi(j)+1,inds}(f(end)+(-totpr+1:0));
end
% figure;subplot(311);imagesc(startAligned);subplot(312);imagesc(endAligned);
startAligned(:,1:find(mean(abs(diff(startAligned'))') > 0,1)-1) = [];
info = info(fi,:);
startAligned(info(:,1) == 3,:) = [];
endAligned(info(:,1) == 3,:) = [];
info(info(:,1) == 3,:) = [];

ilds = [6 4 2 .5];

for i = 1:4
    inds = abs(info(:,2)) <= ilds(i);
    Xs = bsxfun(@rdivide,startAligned(inds,:),std(startAligned(inds,:)));
    Xe = bsxfun(@rdivide,endAligned(inds,:),std(endAligned(inds,:)));
    %Xe = Xe(:,1:size(Xs,2));
    try
    [~,errs(i),~,~,temp] = classify(randn(1,size(Xs,2)),Xs,info(inds,class));%,'diaglinear');
    coeffs(i,:) = temp(1,2).linear;
    catch
        %figure;imagesc(Xs);
    end
    try
    [~,erre(i),~,~,temp] = classify(randn(1,size(Xe,2)),Xe,info(inds,class));%,'diaglinear');
    coeffe(i,:) = temp(1,2).linear;
    catch
        %figure;imagesc(Xe);
    end
    %[~,errb(i),~,~,temp] = classify(randn(1,size(endAligned,2)+size(startAligned,2)),zscore([startAligned(inds,:) endAligned(inds,:)]),info(inds,class));%,'diaglinear');
    %coeffb(i,:) = temp(1,2).linear;
    coeffGs(i,:) = (Xs'*-(info(inds,class)-1.5))'/(Xs'*Xs + lambda*size(Xs,1)*eye(size(Xs,2)));
    yHat = (coeffGs(i,:)*Xs') < 0;
    errGs(i) = sum(yHat'+1 ~=info(inds,class))/sum(inds);
    coeffGe(i,:) = (Xe'*-(info(inds,class)-1.5))'/(Xe'*Xe + lambda*size(Xe,1)*eye(size(Xe,2)));
    yHat = (coeffGe(i,:)*Xe') < 0;
    errGe(i) = sum(yHat'+1 ~=info(inds,class))/sum(inds);
    coeffGb(i,:) = ([Xs Xe]'*-(info(inds,class)-1.5))'/([Xs Xe]'*[Xs Xe] + lambda*size(Xe,1)*eye(size(Xe,2)+size(Xs,2)));
    %yHat = (coeffGe(i,:)*Xe') < 0;
    %errGe(i) = sum(yHat'+1 ~=info(inds,class))/sum(inds);
end

figure;subplot(231);plot(errs);hold all;plot(erre);plot(errGs);plot(errGe);%plot(errb);
subplot(232);plot(coeffs');axis tight;
subplot(233);plot(coeffe');axis tight;
subplot(234);plot(coeffGb');axis tight;
subplot(235);plot(coeffGs');axis tight;
subplot(236);plot(coeffGe');axis tight;

%temp(temp == 0) = nan;

%b = glmfit(temp,info(f,1));
% pupSize = max(1,(squeeze(temp(:,:,5))-40)*2);
% stims = temp(:,:,1);
% inds = [2 3 7 8];
% temp = temp(:,:,inds);
% for i = 1:4
%     t = squeeze(temp(:,:,i));
%     temp(:,:,i) = temp(:,:,i)-nanmean(t(:));
%     temp(:,:,i) = temp(:,:,i)/nanstd(t(:));
% end