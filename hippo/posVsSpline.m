function [cr,cc1,cc2,cm,kerns] = posVsSpline(pos,sp,v,numNeuro)%,vInterp
bounds = [.1 .9];accumbins = 50;
shift = 1;shifts = 0;
if numel(sp) > size(v,1)
    sp = sp(:,1:size(v,1));
end
pos(pos == -1) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i=1:size(pos,2)
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1));
pos = pos(~nanInds,:);v = v(~nanInds,:);sp = sp(:,~nanInds);
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,4);
pos = bsxfun(@minus,pos,mean(pos));
[pos,~,~] = svd(pos(:,1:2),'econ');
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
%v = bsxfun(@rdivide,v,std(v));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
numCross = 3;
%sp = bsxfun(@rdivide,sp,std(sp,0,2));
sp = bsxfun(@times,sp,(v(:,1)./abs(v(:,1))).');
v(:,2) = v(:,1).*conj(v(:,2))./abs(v(:,1));
v(:,1) = [zeros(shift,1); v((shift+1):end,1).*conj(v(1:end-shift,1))./abs(v(1:end-shift,1))];
v = filtLow(v.',1250/32,2).';
sp = filtLow(sp,1250/32,2);
runs = bwlabel(b > 0);
vInterp = zeros(2,2,max(runs),accumbins);
spInterp = zeros(2,size(sp,1),max(runs),accumbins);
bins = (bounds(1))+((1:accumbins)-.5)/accumbins*(diff(bounds));
numX = ceil(sqrt(numNeuro(end)));numY = ceil(numNeuro(end)/numX);
for k = 1:2
    runs = bwlabel(b*((-1)^k)>0);
    runs(vel < .05) = 0;
    [~,spSort(k,:)] = sort(sum(abs(sp(:,runs > 0)),2),'descend');
    for i = 1:max(runs)
        inds = find(runs == i);%inds = min(inds):max(inds);
        %inds(vel(inds,1) < .1) = [];
        for j = 1:2
            vInterp(k,j,i,:) = csaps(pos(inds,1),v(inds,j),1-1e-7,bins);
        end
        for j = 1:numNeuro(end)
            spInterp(k,j,i,:) = csaps(pos(inds,1),sp(spSort(k,j),inds),1-1e-7,bins);
        end
    end
    figure;
    for j = 1:numNeuro(end)
       subplot(numX,numY,j);imagesc(complexIm(squeeze(spInterp(k,j,:,:)),0,1));title(num2str(j));
    end
end
kernFig = figure;
drawnow;h = figure;
shuffle = randperm(max(runs));
cc1 = zeros(numel(numNeuro),2,numCross);cc2 = cc1;cr = cc1;
mc1 = cc1;mc2 = cc2;mr = cr;mm = zeros(2,numCross);cm = mm;
for k = 1:2
    allX = squeeze(spInterp(k,1:size(sp,1),:,:));
    Y = squeeze(vInterp(k,2,shuffle,(shifts+1):(end-shifts)));
    Y1 = squeeze(vInterp(mod(k,2)+1,2,shuffle,shifts+1:end-shifts));
    for n = 1:numel(numNeuro)
        allXC = allX(1:numNeuro(n),:,:);
        [cc1(n,k,:),mc1(n,k,:)] = pipeLine(Y(:).',makeX(allXC,shifts),numCross,10000);
        allXC = allX(1:numNeuro(n),shuffle,:);
        [cc2(n,k,:),mc2(n,k,:)] = pipeLine(Y1(:).',makeX(allXC,shifts),numCross,10000);
        [cr(n,k,:),mr(n,k,:),kern] = pipeLine(Y(:).',makeX(allXC,shifts),numCross,10000);
    end
    [cm(k,:) mm(k,:)] = meanFit(Y,numCross);
%    temp = repmat(mean(Y),[size(Y,1) 1]);
%    cm(k) = corr(Y(:),temp(:));
%    mm(k) = mean(abs(Y(:)-temp(:)).^2)/var(Y(:));%mean(abs(Y(:)).^2);%
    %size(kern)
    kern =  mean(kern,1);kerns(k,:) = squeeze(kern);
    Y = squeeze(vInterp(k,2,:,(shifts+1):(end-shifts)));
    allXC = allX(1:numNeuro(n),:,:);
    allXC = makeX(allXC,shifts);
    subplot(2,4,(k-1)*4+3);image(complexIm(reshape(conj(kern*allXC'),size(Y)),0,.25));
    subplot(2,4,k*4);image(complexIm(Y,0,.25));
    figure(kernFig);
    im = makeKern(kern,numNeuro(end),size(Y,2));
    for i = 1:size(im,3)
        subplot(2,size(im,3),(k-1)*size(im,3)+i);image(complexIm(squeeze(im(:,:,i)),0,1));
    end
    in = squeeze(spInterp(k,1,:,shifts+1:end-shifts));%input('which neuron? ')
    figure(h);subplot(2,4,(k-1)*4+1);imagesc(abs(in));
    subplot(2,4,(k-1)*4+2);image(complexIm(in,0,.25));
end
figure;
subplot(2,4,5);image(complexIm(reshape(conj(kern*allXC'),size(Y)),0,.25));axis off;
subplot(2,4,1);image(complexIm(Y,0,.25));axis off
is = [2 6 7 8 17 18];
pos = [6 7 4 2 3 8];
for i = 1:6
    subplot(2,4,pos(i));image(complexIm(squeeze(spInterp(k,is(i),:,:)),0,1));axis off;
    %subplot(2,4,i+1);imagesc(abs(squeeze(spInterp(k,i,:,:))));axis off;
end
[mean(abs(cr(:)).^2) mean(abs(cm(:)).^2) mean(abs(cc1(:)).^2) mean(abs(cc2(:)).^2)]%;...
[mean(mr(:)) mean(mm(:)) mean(mc1(:)) mean(mc2(:))]

function kerns = makeKern(kern,numNeuro,len)
kern = reshape(kern(2:end),numNeuro,numel(kern(2:end))/numNeuro);
numKnots = 5;
numBlocks  = size(kern,2)/(numKnots+1);
if numKnots >  1
    tau = 1:len;k = 3;
    knots = augknt(linspace(1,max(tau),numKnots),k);
    colmat = spcol(knots,k,brk2knt(tau,1));
else
    colmat = ones(size(X,3),1);
end
for i = 1:numBlocks
kerns(:,:,i) = kern(:,(numKnots+1)*(i-1)+(1:numKnots+1))*colmat.';
end

function [cm mm] = meanFit(y,numCross)
samples = size(y,1);
for i = 1:numCross
    testInds = max(1,ceil([(i-1) i]/numCross*samples));
    testInds = testInds(1):testInds(2);
    trainInds = setdiff(1:samples,testInds);
    if numCross == 1 trainInds = testInds; end
    yTest = y(testInds,:);
    yTrain = y(trainInds,:);
    yEst = repmat(mean(yTrain),[size(yTest,1) 1]);
    cm(i) = corr(yTest(:),yEst(:));
    mm(i) = mean(abs(yTest(:)-yEst(:)).^2)/var(yTest(:));
end

function X1 = makeX(X,shifts)
numKnots = 5;
sx = size(X,2);
for i = 1:shifts
    X = [X; circshift(X(:,1:sx,:),[0 0 -i]) ; circshift(X(:,1:sx,:),[0 0 i]) ];
end
X = X(:,:,(shifts+1):(end-shifts));
if numKnots >  1
    tau = 1:size(X,3);k = 3;
    knots = augknt(linspace(1,max(tau),numKnots),k);
    colmat = spcol(knots,k,brk2knt(tau,1));
else
    colmat = ones(size(X,3),1);
end
X1 = [];
for i = 1:size(colmat,2)
    temp = bsxfun(@times,X,reshape(colmat(:,i),[1 1 size(colmat,1)]));
    X1 = [X1 temp(:,:).'];
end
X1 = [ones(size(X1,1),1) X1 conj(X1) abs(X1)];%real(X1) imag(X1) ];% abs(X1)