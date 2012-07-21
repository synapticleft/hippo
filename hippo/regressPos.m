function regressPos(pos,v,Xf,thresh)

warning off all;
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
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);vel = [0; vel];
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
Xf = [real(Xf);imag(Xf)];Xf = zscore(Xf,0,2);
runs = bwlabel(b > 0);
% vInterp = zeros(2,size(Xf,1),max(runs),accumbins);
% velInterp = zeros(2,max(runs),accumbins);
% velTrace = zeros(2,max(runs),range(timeBins)+1);
% vVel = zeros(2,size(Xf,1),max(runs),range(timeBins)+1);
% bins = (bounds(1))+((1:accumbins)-.5)/accumbins*(diff(bounds));
for k = 1:2
    dat = b*(-1^k)>0 & vel' > thresh;
    rbfPos = gaussPos((1:20)',10,100);%pos(dat,1)
    figure;plot(rbfPos);
%    runs = bwlabel(b*((-1)^k)>0);
%    runs(vel < thresh) = 0;
% for i = 1:max(runs)
%     inds = find(runs == i);inds = min(inds):max(inds);
%     indsa = min(inds)-100:max(inds);
%     start = find(vel(indsa) > thresh,1);
%     start = max(min(indsa)+start-1,-min(timeBins)+1);
%     velTrace(k,i,:) = vel(start+timeBins);
%     inds(vel(inds,1) < .1) = [];
%     for j = 1:size(Xf,1)
%         vVel(k,j,i,:) = Xf(j,start+timeBins);
%         vInterp(k,j,i,:) = csaps(pos(inds,1),Xf(j,inds),1-1e-7,bins);
%     end
%     velInterp(k,i,:) = csaps(pos(inds,1),vel(inds),1-1e-7,bins);
% end
end

function X1 = gaussPos(X,numPoints,numGauss)
X = X - min(X); X = X/max(X);
xGrid = linspace(0,1,numGauss)';
X1 = km_kernel(xGrid,X,'gauss',1/numGauss);