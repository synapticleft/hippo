function [lrV rlV lrAvg rlAvg] = posChunk(pos,sp,v,bounds,gauss)
accumbins = gauss(1);shift = 1;
if numel(sp) > size(v,1)
    sp = sp(:,1:size(v,1));
end
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
nanInds = find(~isnan(pos(:,1)));
pos(:,1) = interp1(nanInds,pos(nanInds,1),1:size(pos,1));
pos(:,2) = interp1(nanInds,pos(nanInds,2),1:size(pos,1));
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);sp = sp(:,~nanInds);
 pos = bsxfun(@minus,pos,mean(pos));%pos = bsxfun(@rdivide,pos,std(pos));
 [a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
%sp = bsxfun(@rdivide,sp,std(sp,0,2));%sp/std(sp(:));
v = bsxfun(@rdivide,v,std(v));
vp11 = v((shift+1):end,1).*conj(v(1:end-shift,1))./abs(v(1:end-shift,1));
vp11 = [zeros(shift,1); vp11];
vp12 = v(:,1).*conj(v(:,2))./abs(v(:,1));
%vp12 = filtLow(vp12,1250/32,2);vp11 = filtLow(vp11,1250/32,2);
%vp12 = vp11;
sp = bsxfun(@times,sp,(v(:,1)./abs(v(:,1))).');
%sp = filtLow(sp,1250/32,2);
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
lrRuns = bwlabel(b>0);rlRuns = bwlabel(b<0);
lrAvg = zeros(size(sp,1),accumbins,max(lrRuns));lrV = zeros(accumbins,max(lrRuns));
rlAvg = zeros(size(sp,1),accumbins,max(rlRuns));rlV = zeros(accumbins,max(rlRuns));
bins = (bounds(1))+((1:accumbins)-.5)/accumbins*(diff(bounds));
runs = bwlabel(b~=0);
cor12 = 0;cor11 = 0;
for i = 1:max(runs)
    inds = find(runs == i);
    cor12 = cor12 + xcov(abs(vp12(inds)),imag(vp11(inds)),500,'coeff');
    %cor11 = cor11 + xcov(vp11,1000,'coeff');
end
figure;plot((cor12)/max(runs));hold all;plot((cor11)/max(runs));
return
for i = 1:max(lrRuns)
    inds = find(lrRuns == i);inds = min(inds):max(inds);
    for j = 1:size(sp,1)
        lrAvg(j,:,i) = csaps(pos(inds,1),sp(j,inds),1-1e-7,bins);
    end
    lrV(:,i) = csaps(pos(inds,1),vp12(inds),1-1e-7,bins);
end
for i = 1:max(rlRuns)
    inds = find(rlRuns == i);inds = min(inds):max(inds);
    for j = 1:size(sp,1)
        rlAvg(j,:,i) = csaps(pos(inds,1),sp(j,inds),1-1e-7,bins);
    end
    rlV(:,i) = csaps(pos(inds,1),vp12(inds),1-1e-7,bins);
end
return
%lrAvg = mean(lrAvg,3);
%lrV = mean(lrV,2);
%rlAvg = mean(rlAvg,3);
%rV = mean(rlV,2);
shifts = 1;
for i = 1:shifts
    lrAvg = [lrAvg; circshift(lrAvg(1:size(sp,1),:,:),[0 -i 0]); circshift(lrAvg(1:size(sp,1),:,:),[0 i 0])];
    rlAvg = [rlAvg; circshift(rlAvg(1:size(sp,1),:,:),[0 -i 0]); circshift(rlAvg(1:size(sp,1),:,:),[0 i 0])];
end
inds = randperm(size(lrAvg,3));
lrV = lrV(:,inds);lrAvg = lrAvg(:,:,inds);
allX = lrAvg(:,:).';%allX = [allX conj(allX) abs(allX)];%ones(size(allX,1),1) 
allX = [real(allX) imag(allX) abs(allX)];
[cc,~,kernlr] = pipeLine(lrV(:).',allX,2,1000);
%figure;subplot(221);imagesc(imag(lrV));
%subplot(222);imagesc(reshape(-imag(mean(kernlr)*allX'),size(lrV)));
temp = mean(reshape(conj(mean(kernlr)*allX'),size(lrV)),2);
temp = temp-mean(temp);temp = temp/std(temp);
lrV = mean(lrV,2);
lrV = lrV-mean(lrV);lrV = lrV/std(lrV);
figure;subplot(211);plot(real(lrV));hold all;plot(imag(lrV));
plot(real(temp));plot(imag(temp));
cc
inds = randperm(size(rlAvg,3));
rlV = rlV(:,inds);rlAvg = rlAvg(:,:,inds);
allX = rlAvg(:,:).';%allX = [allX conj(allX) abs(allX)];%ones(size(allX,1),1) 
allX = [real(allX) imag(allX) abs(allX)];
[cc,~,kernrl] = pipeLine(rlV(:).',allX,2,1000);
temp = mean(reshape(conj(mean(kernrl)*allX'),size(rlV)),2);
temp = temp-mean(temp);temp = temp/std(temp);
rlV = mean(rlV,2);
rlV = rlV - mean(rlV);rlV = rlV/std(rlV);
subplot(212);plot(real(rlV));hold all;plot(imag(rlV));
plot(real(temp));plot(imag(temp));
cc
return
%lrAvg = abs(mean(lrAvg,3));lrAvg = bsxfun(@minus,lrAvg,mean(lrAvg,2));
%lrAvg = bsxfun(@rdivide,lrAvg,std(lrAvg,0,2));
%[u,~,v] = svds(lrAvg,1);
%u = ffa(lrAvg.',1);
%v = pinv(u)*lrAvg;
%figure;subplot(211);plot(real(v));hold all;plot(imag(v));
%rlAvg = abs(mean(rlAvg,3));rlAvg = bsxfun(@minus,rlAvg,mean(rlAvg,2));
%rlAvg = bsxfun(@rdivide,rlAvg,std(rlAvg,0,2));
%[u,~,v] = svds(rlAvg,1);
%u = ffa(rlAvg.',1);
%v = pinv(u)*rlAvg;
%subplot(212);plot(real(v));hold all;plot(imag(v));
%subplot(223);imagesc(imag(rlV));
%subplot(224);imagesc(reshape(-imag(mean(kernrl)*allX'),size(rlV)));
%figure;plot(abs(mean(kernlr)));hold all;plot(abs(mean(kernrl)));
%corr(nanmean(kernlr)',nanmean(kernrl)')
%temp = lrV;lrV = rlV; rlV = temp;
%lrV = rlV;lrAvg = rlAvg;
%lrV = lrV(:,randperm(size(lrV,2)));rlV = rlV(:,randperm(size(rlV,2)));
%lrAvg=abs(lrAvg);lrV = real(lrV);
%%using GLM
trainInds = 1:max(lrRuns);testInds = trainInds;%
trainInds = 1:2:max(lrRuns);%1:floor(max(lrRuns)*.7);
testInds = 2:2:max(lrRuns);%max(trainInds)+1:max(lrRuns);
xTrain = lrAvg(:,:,trainInds);
xTrain = [xTrain(:,:)' xTrain(:,:).' abs(xTrain(:,:)')];
xTest = lrAvg(:,:,testInds);xTest = [xTest(:,:)' xTest(:,:).' abs(xTest(:,:)')];
%[xTrain WM] = whiten(xTrain);xTest = (WM*xTest')';
yTrain = lrV(:,trainInds);yTest = lrV(:,testInds);
a = glmfit(xTrain,yTrain(:));
yHat = glmval(a,xTest,'identity');
corr(yTest(:),yHat)
trainInds = 1:2:max(rlRuns);%1:max(rlRuns);testInds = trainInds;%
testInds = 2:2:max(rlRuns);%1:floor(max(rlRuns)*.7);testInds = max(trainInds)+1:max(rlRuns);
xTrain = rlAvg(:,:,trainInds);
xTrain = [xTrain(:,:)' xTrain(:,:).' abs(xTrain(:,:)')];
xTest = rlAvg(:,:,testInds);xTest = [xTest(:,:)' xTest(:,:).' abs(xTest(:,:)')];
%[xTrain WM] = whiten(xTrain);xTest = (WM*xTest')';
yTrain = rlV(:,trainInds);yTest = rlV(:,testInds);
a1 = glmfit(xTrain,yTrain(:));
figure;plot(abs(a));hold all;plot(abs(a1))
yHat = glmval(a1,xTest,'identity');
yHat1 = glmval(a,xTest,'identity');
%figure;subplot(211);imagesc(real(reshape(lrVHat,size(lrV))));
%subplot(212);imagesc(real(lrV));
%figure;plot(real(lrV(:)));hold all;plot(real(lrVHat));
corr(yTest(:),yHat)
corr(yTest(:),yHat1)
% weighted = bsxfun(@times,rlAvg+randn(size(rlAvg))*.01,gaussWeights);
% [~,indsrl] = sort(mdscale(pdist(abs(weighted)','correlation'),1));%[real(weighted); imag(weighted)]','criterion','sstress'
% weighted = bsxfun(@times,lrAvg+randn(size(lrAvg))*.01,gaussWeights);
% [~,indslr] = sort(mdscale(pdist(abs(weighted)','correlation'),1));%,'criterion','sstress'[real(weighted); imag(weighted)]
% rlAvg = rlAvg(:,indsrl);
% rlV = rlV(:,indsrl);
% rlVa(:,:,1) = (angle(rlV)+pi)/(2*pi);rlVa(:,:,2) = 1;rlVa(:,:,3) = abs(rlAvg)/max(abs(rlAvg(:)));
% rlVa = hsv2rgb(rlVa);
% lrAvg = lrAvg(:,indslr);%lrV = lrV(:,indslr);
% lrV = lrV(:,indslr);
% lrVa(:,:,1) = (angle(lrV)+pi)/(2*pi);lrVa(:,:,2) = 1;lrVa(:,:,3) = abs(lrAvg)/max(abs(lrAvg(:)));
% lrVa = hsv2rgb(lrVa);

function [X WM] = whiten(X)
A = X'*X/size(X,1);
[V,D] = eig(A);d = diag(D);
D1 = diag(sqrt(1./d));
WM = D1*V';
X = (WM*X')';
DWM = V*D1^(-1);