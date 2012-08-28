function regressPos(pos,v,Xf,thresh)

warning off all;
s = size(Xf,1);
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
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);%sp = sp(:,~nanInds);
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);vel = [0; vel]/max(vel);
nanInds = find(~isnan(vel));
vel = interp1(nanInds,vel(nanInds),1:numel(vel));
pos = bsxfun(@minus,pos,mean(pos));%pos = bsxfun(@rdivide,pos,std(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
pos(:,1) = pos(:,1)-min(pos(:,1));pos(:,1) = pos(:,1)/max(pos(:,1));
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
%v(:,2) = v(:,2).*conj(v(:,1))./abs(v(:,1));
offSet = 1;
%v(:,1) = [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))];
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).'); ...
    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xf = [real(Xf);imag(Xf)];
%Xf = gaussPos(pos(:,1),10);
Xf = zscore(Xf,0,2);
%figure;plot(Xf');return
runs = bwlabel(b > 0);Xf1 = [];
bases = 10;
kern = zeros(2,bases,size(Xf,1));
Xf1 = zeros(2*bases,size(Xf,2));
link = 'identity';%'logit';
numCross = 3;c = [];
for k = 1:2
    %dat = b*((-1)^k)>0 & vel > thresh;
    data = bwlabel(b*((-1)^k)>0);
    vec = 1:max(data);randTrials = ceil(rand(size(vec))*numCross);
    %rbfPos = gaussPos(pos(dat,1),bases);
    %rbfTrial = gaussPos(data(dat)',bases);
    randTrials = randperm(max(data));
    for i = 1:numCross
        testInds = vec(randTrials == i);%max(1,ceil([(i-1) i]/numCross*samples));
        %testInds = testInds(1):testInds(2);
        trainInds = vec(randTrials ~= i);%setdiff(1:samples,testInds);
        dat = b*((-1)^k)>0 & vel > thresh & ismember(data,trainInds);
        rbfPos = gaussPos(pos(dat,1),bases);%
        %XrbfPos = gaussPos(data(dat)',bases);
        tempKern(i,:,:) = ((Xf(:,dat)*Xf(:,dat)'+1000)\Xf(:,dat)*rbfPos')';
        dat = b*((-1)^k)>0 & vel > thresh & ismember(data,testInds);
        rbfPos = gaussPos(pos(dat,1),bases);
        %rbfPos = gaussPos(data(dat)',bases);
        temp = squeeze(tempKern(i,:,:))*Xf(:,dat);
        c = [c corr(rbfPos(:),temp(:))];
    end
    %pos(dat,1)
%    for i = 1:size(rbfPos,1)
%        kern(k,i,:) = glmfit(Xf(:,dat).',rbfPos(i,:)','normal',link);
%        Xf1((k-1)*bases+i,:) = glmval(squeeze(kern(k,i,:)),Xf',link)';
%    end        
    kern(k,:,:) = mean(tempKern);
    Xf1((k-1)*bases+(1:bases),:) = squeeze(kern(k,:,:))*Xf;
end
Xf = Xf1;clear Xf1;
% figure;subplot(211);plot(Xf(1:10,:)');subplot(212);plot(Xf(11:20,:)');
% sPlot(Xf);return
vInterp = zeros(2,size(Xf,1),max(runs),accumbins);
for k = 1:2
   runs = bwlabel(b*((-1)^k)>0);
   runs(vel < thresh) = 0;
   for i = 1:max(runs)
       inds = find(runs == i);inds = min(inds):max(inds);
%       indsa = min(inds)-100:max(inds);
%       start = find(vel(indsa) > thresh,1);
%       start = max(min(indsa)+start-1,-min(timeBins)+1);
%       velTrace(k,i,:) = vel(start+timeBins);
%       inds(vel(inds) < .1) = [];
       for j = 1:size(Xf,1)
%           vVel(k,j,i,:) = Xf(j,start+timeBins);
           %vInterp(k,j,i,:) = csaps(pos(inds,1),Xf(j,inds),1-1e-7,bins);
           vInterp(k,j,i,:) = accumarray(max(1,min(accumbins,floor((pos(inds,1)-bounds(1))...
               /(bounds(2)-bounds(1))*accumbins)+1)),Xf(j,inds),[accumbins 1],@mean);
       end
%       velInterp(k,i,:) = csaps(pos(inds,1),vel(inds),1-1e-7,bins);
   end
end
t = [squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))];
xdim = ceil(sqrt(size(t,1)));ydim = ceil(size(t,1)/xdim);
figure;for i = 1:size(t,1)
subplot(xdim,ydim,i);imagesc(imfilter(reshape(t(i,:),[max(runs) 2*accumbins]),fspecial('gaussian',5,1)));axis off;
end
superImp(reshape([squeeze(vInterp(1,:,:)) squeeze(vInterp(2,:,:))],...
   [size(vInterp,2) size(vInterp,3) 2*size(vInterp,4)]),1:size(vInterp,2),1);
xdim = ceil(sqrt(2*size(kern,2)));ydim = ceil(2*size(kern,2)/xdim);
figure;for k = 1:2
    k1 = squeeze(kern(k,:,:));
    %k1 = pinv(k1)';
    for i = 1:size(kern,2)
        subplot(xdim,ydim,(k-1)*size(kern,2)+i);imagesc(complexIm(reshape(...
            squeeze(complex(k1(i,1:s),k1(i,(s+2):(end-1)))),[8 s/8]),0,1));
    end
end

function X1 = gaussPos(X,numPoints)
X = X - min(X); X = X/max(X);
xGrid = linspace(0,1,numPoints)';
X1 = km_kernel(xGrid,X,'gauss',(1/(2*numPoints)));