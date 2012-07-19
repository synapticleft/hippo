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
vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);
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