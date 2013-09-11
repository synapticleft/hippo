function sessionsXcor(onlyCA1)
%%plot decay of coding similarity across sessions, compare:
% 1) spikes vs. spf vs. LFP
% 2) MSE vs. xcov vs. xcorr
% 3) original vs. whitened
% 4) within sessions vs. across sessions  %%maybe not across sessions since
% track was rotated

ratID = 'ec013';
dataDir = '/media/Expansion Drive/gaMatlab/';
load('/media/Expansion Drive/redwood/KenjiData.mat','Beh');
d = dir([dataDir ratID '*.*']);
%inds = [33:61 63:64];
if onlyCA1
    inds = 33:64;%[1:61 63:64];
    shanks = 5:8;
else
    inds = 1:64;
    shanks = 1:8;
end
nbins = 100;
fLast = 1;
for i = 1:numel(d)
    load([dataDir d(i).name],'pos','Xf','spf','v','cellInfo');
    Xf = Xf(inds,:);
    Xf = bsxfun(@minus,Xf,mean(Xf,2));
    Xf = bsxfun(@rdivide,Xf,std(Xf,0,2));
    pos = pos(1:size(Xf,2),:);%
    [~,pos,thresh,run] = fixPos(pos);
    pos(pos > 1) = 3-pos(pos > 1);
    pos = pos*pi;
    %if i == 1
    [u,~] = eig(Xf(:,thresh)*Xf(:,thresh)');
    u = bsxfun(@times,u,exp(-1i*angle(mean(u))));
    %end
    v = (u(:,1)\Xf)';
    spf = spf(:,1:size(Xf,2));
    spInds = ismember(cellInfo.shank,shanks) && cellInfo.type;
    spf1 = zeros(cellInfo.popSize,size(spf,2));
    spf1(cellInfo.ID(spInds),:) = spf(spInds,:);
    spf = spf1; clear spf1;
    Xf = bsxfun(@times,Xf,exp(1i*angle(v.')));%spf
    spf = bsxfun(@times,spf,exp(1i*angle(v.')));
    [~,yRbf] = get1Dbasis('vonmises',nbases,pos,70);
    Xf = Xf.';
    %Xf = zscore(Xf);
    f = find(strcmp(Beh(:,4),d(i).name(1:end-4)));
    rn = randperm(sum(thresh));
    trInds = thresh;%(rn(1:floor(sum(thresh)/2)));
    teInds = thresh;%(rn(ceil(sum(thresh)/2):end));
    if ~strcmp(Beh(f,2),Beh(fLast,2)) || ~generalize%size(W,1) ~= size(Xf,2)
        %W = Xf(thresh,:)\yRbf(thresh,:);
        W = (Xf(trInds,:)'*Xf(trInds,:) + numel(trInds)*eye(size(Xf,2)))\(Xf(trInds,:)'*yRbf(trInds,:));
    end
    fLast = f;
%    for j = 1:size(Xf,2)
%        Xfd(:,j) = decimate(Xf(:,j),4);
%    end
%    Xf = decimate(Xf.',4);
    yhat = Xf*W;
    for j = 1:nbases
        yhPos(j,:) = accumarray([min(nbins,ceil(pos(teInds)/2/pi*nbins+eps))],yhat(teInds,j),[nbins 1],@mean);
    end
    figure(1);subplot(5,5,i);imagesc(abs(yhPos));axis off;title(Beh(f,2));
    if useSpikes
    xPos = zeros(size(Xf,2),max(run),nbins);
    for j = 1:size(Xf,2)
       xPos(j,:,:) = accumarray([run(thresh)' min(nbins,ceil(pos(thresh)/2/pi*nbins+eps))],Xf(thresh,j),[max(run) nbins],@mean);
       figure;showGrid(xPos,[],.5);title(Beh(f,2));
    end
    else
        xPos = zeros(size(Xf,2),nbins);
        for j = 1:size(Xf,2)
            xPos(j,:) = accumarray([min(nbins,ceil(pos(thresh)/2/pi*nbins+eps))],Xf(thresh,j),[nbins 1],@mean);
        end
        figure(2);subplot(5,5,i);imagesc(complexIm(xPos,0,2,4));axis off;title(Beh(f,2));
        xPos = bsxfun(@rdivide,xPos,sqrt(sum(xPos.*conj(xPos))));
        xPos = bsxfun(@minus,xPos,mean(xPos.').');
        figure(3);subplot(5,5,i);imagesc(complexIm(xPos,0,1,1));axis off;title(Beh(f,2));
    end
    %subplot(5,5,i);imagesc(complexIm(reshape(u(:,2),[8 4]),0,2,5));title(cellInfo.popSize);drawnow;
    %figure(2);subplot(3,3,i);imagesc(abs(xPos));drawnow;
%    [~,m] = max(abs(yhat)');
%    figure;imagesc(abs(yhat(thresh,:))');
%    hold all;scatter(1:sum(thresh),m(thresh),'r','filled');
end