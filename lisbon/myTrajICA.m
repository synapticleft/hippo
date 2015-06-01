function [wx Xor trialNum stimOr X] = myTrajICA(fn,ILD)

inds = [15 16];
if ~exist('ILD','var')
    ILD = [.5 2 4];
end
%[linspace(0,1,64);linspace(0,1,64);linspace(1,0,64)]';
%cols = colormap('jet');
%d = preProcessRoberto(fn,[inds 6],0,0,[],[-ILD ILD]);

file = dir('*.mat');
load(file(fn).name,'data','center_ILD');
choice = [data{2:end,7}];
f = choice ~= 3;

scales = 2.^(-1:.2:2);
X = [];Xor = [];stim = [];stimOr = [];
version = 2;
for i = 2:size(data,1)%1:size(d,1)
    if f(i-1)
    %f = squeeze(d(i,:,end)) ~= 0;
    %f(end-10:end) = 0;
    trialData = reshape([data{i,inds}],[numel(data{i,6}) numel(inds)])';
    traj = trialData(1,:) + 1i*trialData(2,:);
    %trialData = trialData(:,~isnan(data{i,6}));
    %traj = squeeze(d(i,:,1) + 1i*d(i,:,2));%
    if version == 1
        trajcwt = [cwt(real(diff(traj,[],2)),scales,'cgau1'); cwt(imag(diff(traj,[],2)),scales,'cgau1')];
        trajcwt = [zeros(size(trajcwt,1),1) trajcwt];
    elseif version == 2
            trajcwt = [cwt(real(traj),scales,'cgau1'); cwt(imag(traj),scales,'cgau1')];
    else
            trajcwt = cwt(traj,scales,'cgau1');
    end
    trajcwt = trajcwt(:,~isnan(data{i,6}));
    %trajcwt = [cwt(squeeze(d(i,:,1)),scales,'cgau1');cwt(squeeze(d(i,:,2)),scales,'cgau1')];
    %trajcwt = log(max(eps,abs(cwt(squeeze(d(i,:,1)),scales,'cgau1')))) + 1i*log(max(eps,abs(cwt(squeeze(d(i,:,2)),scales,'cgau1'))));
    stimDat = data{i,6} - center_ILD;
    stimcwt = cwt(stimDat(~isnan(data{i,6})),scales,'cgau1');
    X = [X trajcwt];%(:,f)];
    Xor = [Xor traj(:,~isnan(data{i,6}))];%(f)];
    stim = [stim stimcwt];
    stimOr = [stimOr stimDat(~isnan(data{i,6}))];
    %dx = squeeze(d(i,f,1));
    %dy = squeeze(d(i,f,2));
    %dx = filtfilt(gausswin(8),sum(gausswin(8)),dx);
    %dy = filtfilt(gausswin(8),sum(gausswin(8)),dy);
    %for j = ord:numel(f)-precludeEnd
    %    X =[X [dx(j+(-ord+1:0))'; dy(j+(-ord+1:0))']];
    %    y =[y dx(j+(predictFwd))' > dx(j)];
    %end
    trialNum(i) = size(X,2);
    end
%    input('');
end
%f = find(f);
%X = bsxfun(@rdivide,X,std(X,0,2));
isReal = 0;
if isReal
[wx A W ] = fastica([real(X);imag(X)],'approach','symm');
else
[A W] = cfastica(X);%,'mle_circ');% Z
%[U,S,V] = svd(X,'econ');
%A = U*S;W = pinv(A);
wx = W*X;
%wx = bsxfun(@rdivide,wx,std(wx,0,2));
for i = 1:size(wx,1)
    [temp,~] = svds([real(wx(i,:));imag(wx(i,:))],1);
    ang(i) = exp(1i*atan(temp(2)/temp(1)));%angle(mean(sign(real(wx)).*wx,2)).';
end
temp = [];
%ang = angle(mean(A(1:end/2,:)));
A = bsxfun(@times,A,ang);W = bsxfun(@times,W,conj(ang).');wx = bsxfun(@times,wx,conj(ang).');
end
figure;scatter(real(A(:)),imag(A(:)),'filled');drawnow;
[~,s] = sort(sum(A.*conj(A)),'descend');
A = A(:,s); W = W(s,:);wx = wx(s,:);
figure;imagesc(complexIm(A));
%figure;imagesc(complexIm(A));
trialNum = trialNum(trialNum ~= 0);
trialNum = [0 trialNum];
%return
figure;
xs{1} = -3:.1:3;xs{2} = xs{1};for i = 1:min(16,size(wx,1))
subplot(4,4,i);imagesc(log(hist3([real(wx(i,:)); imag(wx(i,:))]',xs)));
end
%wx = real(wx);
%return
figure;for i = 1:min(25,size(wx,1))
subplot(5,5,i);
for j = 1:size(wx,1)
temp(j,:) = ((xcov(wx(i,:),wx(j,:),100,'coeff')));
end
imagesc(complexIm(temp));
end
drawnow;

clear temp;

figure;for i = 1:min(16,size(wx,1))
subplot(4,4,i);
for j = 1:size(stim,1)
temp(j,:) = ((xcov(wx(i,:),stim(j,:),100,'coeff')));
end
imagesc(complexIm(temp));
end
%return
figure;
cols = colormap('jet');
%for i = 1:numel(trialNum)-1
    inds = trialNum(30):trialNum(50);%trialNum(i)+1:trialNum(i+1);
    for j = 1:min(25,size(wx,1))
        %acc1 = abs(W(j,:)*X(:,inds));acc1 = [max(1,min(ceil(acc1/max(acc1)*64),64))];
        act = wx(j,inds);
         indPlot = abs(act) > 1;
         act = act(indPlot);%act = wx(j,inds(indPlot));
        actI = max(1,min(64,32+round(10*real(act))));
        %cols = squeeze(complexIm(act,0,.5));
        subplot(5,5,j);scatter(real(Xor(inds(indPlot))),imag(Xor(inds(indPlot))),abs(act)*10,cols(actI,:),'filled');%hold all;
        axis tight
        drawnow;
    end
    %input('');
%end
%     acc = X([ord end],inds); %- X([ord end]-1,inds);%*2 + X([ord end] -2,inds);
%     fs = scal2frq(2.^(0:.1:5),'cgau1',1/80);
%     tempx = cwt(acc(1,:),2.^(0:.1:5),'cgau1');
%     tempy = cwt(acc(2,:),2.^(0:.1:5),'cgau1');
%     tempz = cwt(acc(1,:) + 1i*acc(2,:),2.^(0:.1:5),'cgau1');
%     temp = sqrt(abs(tempx).^2 + abs(tempy).^2);%temp = temp(:,2:end-1);
%     tempa = abs(tempx + tempy);
%     subplot(323);imagesc([1 size(acc,2)],fs,complexIm(bsxfun(@rdivide,tempx,max(abs(tempx),[],2)),0,1));
%     subplot(325);imagesc([1 size(acc,2)],fs,complexIm(bsxfun(@rdivide,tempz,max(abs(tempy),[],2)),0,1));
%     subplot(324);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;
%     subplot(326);plot(sqrt(bsxfun(@rdivide,tempa',max(abs(tempz')))));axis tight;
%     %subplot(323);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;%imagesc(1:size(temp,2),fs,complexIm(bsxfun(@rdivide,temp,max(abs(temp),[],2))));%plot(acc);
%     %subplot(324);imagesc(1:size(temp,1),fs,sqrt(bsxfun(@rdivide,temp',max(abs(temp'))))');
%     
% %      fs = scal2frq(2.^(-1:.1:5),'cgau2',1/80);
% %     temp = abs(cwt(acc(1,:),2.^(-1:.1:5),'cgau2')).^2 + abs(cwt(acc(2,:),2.^(-1:.1:5),'cgau2')).^2;%temp = temp(:,2:end-1);
% %     subplot(325);plot(sqrt(bsxfun(@rdivide,temp',max(abs(temp')))));axis tight;%imagesc(1:size(temp,2),fs,complexIm(bsxfun(@rdivide,temp,max(abs(temp),[],2))));
% %     subplot(326);imagesc(1:size(temp,1),fs,sqrt(bsxfun(@rdivide,temp',max(abs(temp'))))');
% %     
%     acc(1,:) = filtfilt(gausswin(8),sum(gausswin(8)),acc(1,:));
%     acc(2,:) = filtfilt(gausswin(8),sum(gausswin(8)),acc(2,:));
%     acc1 = sqrt(sum(diff(acc,[],2).^2));acc1 = [1 max(1,min(ceil(acc1/max(acc1)*64),64))];
%     subplot(321);scatter(1:numel(inds),X(ord,inds),[],cols(acc1,:),'filled');hold on;
%     scatter(1:numel(inds),X(end,inds),[],cols(acc1,:),'filled');hold off;axis tight;%ceil(posterior(inds)*64)
%     subplot(322);scatter(1:numel(inds),X(ord,inds),[],cols(acc1,:),'filled');hold on;
%     scatter(1:numel(inds),X(end,inds),[],cols(acc1,:),'filled');hold off;axis tight;%ceil(posterior(inds)*64)
%         input('');
% end    