function tempx = myAffine(fn,ILD,ord)

inds = [15 16];
predictFwd = 2;
precludeEnd = predictFwd;
if ~exist('ILD','var')
    ILD = [.5 2 4];
end
cols = [linspace(0,1,64);zeros(1,64);linspace(1,0,64)]';
%cols = colormap('jet');
[d,~,~,center_ILD] = preProcessRoberto(fn,[inds 6],0,0,[],[-ILD ILD]);

xs{2} = -1:.05:2;xs{1} = -3:.05:3;

h = hist3([makeFlat(d(:,:,1)) makeFlat(d(:,:,2))],xs);
h = imfilter(h,fspecial('gaussian',10,3));
ma = max(max(h(xs{1} < 0,xs{2} > 0)));
mb = max(max(h(xs{1} > 0,xs{2} > 0)));
[xa,ya] = find(h == ma);
[xb,yb] = find(h == mb);
ma = [xs{1}(xa) xs{2}(ya)];
mb = [xs{1}(xb) xs{2}(yb)];

X = [];y = [];
for i = 1:size(d,1)
    f = find(squeeze(d(i,:,end)) ~= 0);
    dx = squeeze(d(i,f,1));
    dy = squeeze(d(i,f,2));
    %dx = filtfilt(gausswin(8),sum(gausswin(8)),dx);
    %dy = filtfilt(gausswin(8),sum(gausswin(8)),dy);
    for j = ord:numel(f)-precludeEnd
        X =[X [dx(j+(-ord+1:0))'; dy(j+(-ord+1:0))']];
        y =[y dx(j+(predictFwd))' > dx(j)];
%        transformPointsForward(temp,[dx(j+1+(-ord:0)); dy(j+1+(-ord:0))]');
    end
    trialNum(i) = numel(y);
%    input('');
end
trialNum = [0 trialNum];
% for i = 1:size(y,1)
%     [class(i,:),err(i),p,~,c] = classify(X',X',y(i,:)+1);
%     posterior(i,:) = p(:,2);
%     coeff(i,:) = [c(2,1).linear' c(2,1).const];
% end

figure;for i = 1:size(d,1)
    inds = trialNum(i)+1:trialNum(i+1);
    acc = X([ord end],inds); %- X([ord end]-1,inds);%*2 + X([ord end] -2,inds);
    fs = scal2frq(2.^(0:.1:5),'cgau1',1/80);
    tempx = cwt(acc(1,:),2.^(0:.1:5),'cgau1');
    tempy = cwt(acc(2,:),2.^(0:.1:5),'cgau1');
    tempz = cwt(acc(1,:) + 1i*acc(2,:),2.^(0:.1:5),'cgau1');
    temp = sqrt(abs(tempx).^2 + abs(tempy).^2);%temp = temp(:,2:end-1);
    tempa = abs(tempx + tempy);
    subplot(323);imagesc([1 size(acc,2)],fs,complexIm(bsxfun(@rdivide,tempx,max(abs(tempx),[],2)),0,1));
    subplot(325);imagesc([1 size(acc,2)],fs,complexIm(bsxfun(@rdivide,tempz,max(abs(tempy),[],2)),0,1));
    subplot(324);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;
    subplot(326);plot(sqrt(bsxfun(@rdivide,tempa',max(abs(tempz')))));axis tight;
    %subplot(323);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;%imagesc(1:size(temp,2),fs,complexIm(bsxfun(@rdivide,temp,max(abs(temp),[],2))));%plot(acc);
    %subplot(324);imagesc(1:size(temp,1),fs,sqrt(bsxfun(@rdivide,temp',max(abs(temp'))))');
    
%      fs = scal2frq(2.^(-1:.1:5),'cgau2',1/80);
%     temp = abs(cwt(acc(1,:),2.^(-1:.1:5),'cgau2')).^2 + abs(cwt(acc(2,:),2.^(-1:.1:5),'cgau2')).^2;%temp = temp(:,2:end-1);
%     subplot(325);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;%imagesc(1:size(temp,2),fs,complexIm(bsxfun(@rdivide,temp,max(abs(temp),[],2))));
%     subplot(326);imagesc(1:size(temp,1),fs,sqrt(bsxfun(@rdivide,temp',max(abs(temp'))))');
%     
    acc(1,:) = filtfilt(gausswin(8),sum(gausswin(8)),acc(1,:));
    acc(2,:) = filtfilt(gausswin(8),sum(gausswin(8)),acc(2,:));
    acc1 = sqrt(sum(diff(acc,[],2).^2));acc1 = [1 max(1,min(ceil(acc1/max(acc1)*64),64))];
    subplot(321);scatter(1:numel(inds),X(ord,inds),[],cols(acc1,:),'filled');hold on;
    scatter(1:numel(inds),X(end,inds),[],cols(acc1,:),'filled');hold off;axis tight;%ceil(posterior(inds)*64)
    subplot(322);scatter(1:numel(inds),X(ord,inds),[],cols(acc1,:),'filled');hold on;
    scatter(1:numel(inds),X(end,inds),[],cols(acc1,:),'filled');hold off;axis tight;%ceil(posterior(inds)*64)
        input('');
end    
%figure;imagesc(corr(posterior'));
return
%imagesc(y/X);
%figure;plot(y/X);
%figure;plot((y/X)');

function trajxhat = myMin(trajx) %,costHat [mse,cost]
    trajxhat = minimumJerk(trajx(3),trajx(3)-trajx(2),trajx(3)-2*trajx(2)+trajx(1),trajx(end),trajx(end)-trajx(end-1),trajx(end)-2*trajx(end-1)+trajx(end-2),numel(trajx)-2);
%     mse(1) = sum((trajx(3:end) - trajxhat).^2);
%     cost(1) = numel(trajx);%sum(trajx(3:end).^2);
%     for i = 1:3
%         mse(i+1) = sum((diff(trajx(3:end),i) - diff(trajxhat,i)).^2);
%         cost(i+1) = sum(diff(trajx(3:end),i).^2);
%         %costHat(i) = sum(diff(trajxhat,i).^2);
%     end