function decodeAll(ratID,useSpikes,generalize,rev,dec,subSet)
%%check how OLE decoder trained on one session generalizes to other sessions

dataDir = '/media/Expansion Drive/gaMatlab/';
load('/media/Expansion Drive/redwood/KenjiData.mat','Beh');
d = dir([dataDir ratID '*.*']);
%inds = [33:61 63:64];
inds = 1:64;%33:64;%[1:61 63:64];
%inds = 1:32;
nbins = 100;
basis.n = 50;
basis.s = 70;
yhPos = zeros(basis.n,nbins);
fLast = 1;
i = 0;
figure;

% Positions to use for decoding...
pvec = linspace(0,2*pi,nbins); % track length is mapped 0-2*pi for von-mises/fourier bases
[~,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);

finds = 1:numel(d);
if rev
    finds = finds(end:-1:1);
end

for jj = finds%numel(d):-1:
    if ismember(d(jj).name,subSet)
        i = i + 1;
    load([dataDir d(jj).name],'pos','Xf','spf','v','cellInfo');
    Xf = Xf(inds,:);
    pos = pos(1:size(Xf,2),:);%
    vel = angVel(pos);vel = filtLow(vel(:,1),1250/32,1);vel = vel/max(vel);
    [~,pos,thresh,~] = fixPos(pos);
    pos(pos > 1) = 3-pos(pos > 1);
    pos = pos*pi;
    %if i == 1
    [u,~] = eig(Xf(:,thresh)*Xf(:,thresh)');
    u = bsxfun(@times,u,exp(-1i*angle(mean(u))));
    %end
    v = (u(:,1)\Xf)';
    if useSpikes
        spf = spf(:,1:size(Xf,2));
        Xf = zeros(cellInfo.popSize,size(spf,2));
        Xf(cellInfo.ID,:) = spf;
    end
    Xf = bsxfun(@times,Xf,exp(1i*angle(v.')));
    if dec > 1
    for j = 1:size(Xf,1)
        Xfd(j,:) = decimate(Xf(j,:),dec);
    end
    Xf = Xfd;clear Xfd;
    vel = decimate(vel,dec);
    pos = angle(-decimate(exp(1i*pos),dec))+pi;
    thresh = logical(round(decimate(double(thresh),dec)));
    end
    rn = randperm(sum(thresh));
    trInds = thresh;%(rn(1:floor(sum(thresh)/2)));
    teInds = thresh;%(rn(ceil(sum(thresh)/2):end));
    %Xf(:,trInds) = zscore(Xf(:,trInds),0,2);
    [~,yRbf] = get1Dbasis('vonmises',basis.n,pos,basis.s);
    Xf = [real(Xf);imag(Xf)];
    Xf = Xf.';
    f = find(strcmp(Beh(:,4),d(jj).name(1:end-4)));
    if i ==1 || ~generalize %~strcmp(Beh(f,2),Beh(fLast,2)) %size(W,1) ~= size(Xf,2)
        %W = Xf(trInds,:)\yRbf(trInds,:);
        W = (Xf(trInds,:)'*Xf(trInds,:) + numel(trInds)*eye(size(Xf,2)))\(Xf(trInds,:)'*yRbf(trInds,:));
        %W = (Xf'*Xf + size(Xf,1)*eye(size(Xf,2)))\(Xf'*yRbf);
    end
    fLast = f;
    yhat = Xf*(W*dbasis');
    yhat = bsxfun(@rdivide,yhat,sqrt(mean(abs(yhat(teInds,:))).^2));
    %yhat = exp(yhat);
    %yhat = bsxfun(@rdivide,yhat,sum(yhat,2));
    [~,maxpost]=max(abs(yhat)');
    subplot(2,3,3*(i-1)+3);imagesc((hist3([pos(teInds) maxpost(teInds)'],[30 30])));xlabel('predicted');ylabel('actual');
    for j = 1:size(yhat,2)
        yhPos(j,:) = accumarray([min(nbins,ceil(pos(teInds)/2/pi*nbins+eps))],yhat(teInds,j),[nbins 1],@mean);
    end
    someHists(2*i-1,:) = hist(pos(trInds)/2/pi,linspace(0,1,nbins));
    someHists(2*i,:) = accumarray([min(nbins,ceil(pos(teInds)/2/pi*nbins+eps))],vel(trInds),[nbins 1],@mean);
    subplot(2,3,3*(i-1)+1);imagesc(abs(yhPos));axis off;title(Beh(f,2));
    subplot(2,3,3*(i-1)+2);imagesc(abs(bsxfun(@rdivide,yhPos,sqrt(sum(abs(yhPos).^2,2)))));
    err = circ_dist(pos(teInds)*pi,maxpost(teInds)'/size(yhat,2)*2*pi);
    median(abs(err))
     end
end
%someHists = bsxfun(@rdivide,someHists,max(someHists')');
%figure;plot(someHists');