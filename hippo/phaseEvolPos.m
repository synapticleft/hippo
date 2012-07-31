function phaseEvolPos(pos,v,Xf,u,thresh,mName)

bounds = [.1 .9];
accumbins = 200;
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
%v(:,1) = [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))];
offSet = 0;
Xf = bsxfun(@times,Xf,exp(1i*angle(mean(u(:,1)*[zeros(offSet,1) v(1:end-offSet,1).']))));
Xf = complex(zscore(real(Xf),0,2),zscore(imag(Xf),0,2));
%Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1))).'); ...
%    [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];

c = getCol(size(Xf,1));
c1 = getCol(accumbins);
vInterp = zeros(2,size(Xf,1),accumbins);
for k = 1:2
    inds = b*((-1)^k)>0 & vel > thresh;
    %figure;plot(linspace(0,1,20),hist(pos(inds,1),linspace(0,1,20)));hold all;
    inds = bwmorph(inds,'dilate',20);
    %plot(linspace(0,1,20),hist(pos(inds,1),linspace(0,1,20)));return
    bounds = [min(pos(inds)) max(pos(inds))];
    bb = linspace(bounds(1),bounds(2),accumbins);
    for j = 1:size(Xf,1)
        vInterp(k,j,:) = accumarray(max(1,min(accumbins,floor((pos(inds,1)-bounds(1))...
            /(bounds(2)-bounds(1))*accumbins)+1)),Xf(j,inds),[accumbins 1],@mean);
    end
    vInterp(k,:,:) = filtLow(squeeze(vInterp(k,:,:)),4,1);
    figure;
    for i = 1:accumbins
        plot(squeeze(real(vInterp(k,:,i))),squeeze(imag(vInterp(k,:,i))),'color',c1(i,:));hold on;
        s = scatter(squeeze(real(vInterp(k,:,i))),squeeze(imag(vInterp(k,:,i))),60,c,'o','filled');%hold off;
        title(round(100*bb(i)));
        set(gca,'xticklabels',[],'yticklabels',[],'fontsize',20);axis tight;
        drawnow;
        m(i) = getframe(gcf);
        delete(s);
    end
    if exist('mName','var')
        movie2avi(m,[mName num2str(k) '.avi'],'fps',10);
    end
    figure;image(bb,1:size(vInterp,2),complexIm(squeeze(vInterp(k,:,:)),0,1));
end

function c = getCol(n)
c = repmat(linspace(0,1,n)',[1 3]);
for i = 1:3
c(:,i) = max(0,min((1-abs(c(:,i)-i/4)*2.5)*1.5,1));
end
c = fliplr(c);