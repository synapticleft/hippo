function [wx Xor trialNum stimOr X] = myTrajICA1(fn,ILD)

inds = [15];%[15 16];%15 16
stimInd = 6;%6
if ~exist('ILD','var')
    ILD = [.5 2 4 6];
end
%d = preProcessRoberto(fn,[inds 6],0,0,[],[-ILD ILD]);

%file = dir('*.mat');
load(fn,'data','center_ILD');
choice = [data{2:end,7}];
f = choice ~= 3;% & (1:1339 > 315);

scales = 2.^(-1:.2:2);
X = [];Xor = [];stim = [];stimOr = [];
version = 2;
for i = 2:size(data,1)%1:size(d,1)
    if f(i-1)
    %f = squeeze(d(i,:,end)) ~= 0;
    %f(end-10:end) = 0;
    trialData = reshape([data{i,inds}],[numel(data{i,stimInd}) numel(inds)])';
    traj = trialData(1,:);% + 1i*trialData(2,:);
    %trialData = trialData(:,~isnan(data{i,6}));
    %traj = squeeze(d(i,:,1) + 1i*d(i,:,2));%
    if version == 1
        trajcwt = [cwt(real(diff(traj,[],2)),scales,'cgau1'); cwt(imag(diff(traj,[],2)),scales,'cgau1')];
        trajcwt = [zeros(size(trajcwt,1),1) trajcwt];
    elseif version == 2
            trajcwt = [cwt(real(traj),scales,'cgau1')];% cwt(imag(traj),scales,'cgau1')];%
    else
            trajcwt = cwt(traj,scales,'cgau1');
    end
    trajcwt = trajcwt(:,~isnan(data{i,stimInd}));
    %trajcwt = [cwt(squeeze(d(i,:,1)),scales,'cgau1');cwt(squeeze(d(i,:,2)),scales,'cgau1')];
    %trajcwt = log(max(eps,abs(cwt(squeeze(d(i,:,1)),scales,'cgau1')))) + 1i*log(max(eps,abs(cwt(squeeze(d(i,:,2)),scales,'cgau1'))));
    stimDat = data{i,stimInd} - center_ILD;
    stimcwt = cwt(stimDat(~isnan(data{i,stimInd})),scales,'cgau1');
    X = [X trajcwt];%(:,f)];
    Xor = [Xor traj(:,~isnan(data{i,stimInd}))];%(f)];
    stim = [stim stimcwt];
    stimOr = [stimOr stimDat(~isnan(data{i,stimInd}))];
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
[A W] = cfastica(X);
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

% figure;for i = 1:min(16,size(wx,1))
% subplot(4,4,i);
% for j = 1:size(stim,1)
% temp(j,:) = ((xcov(wx(i,:),stim(j,:),100,'coeff')));
% end
% imagesc(complexIm(temp));
% end
% figure;
% cols = colormap('jet');
% %for i = 1:numel(trialNum)-1
%     inds = trialNum(30):trialNum(50);%trialNum(i)+1:trialNum(i+1);
%     for j = 1:min(25,size(wx,1))
%         %acc1 = abs(W(j,:)*X(:,inds));acc1 = [max(1,min(ceil(acc1/max(acc1)*64),64))];
%         act = wx(j,inds);
%          indPlot = abs(act) > 1;
%          act = act(indPlot);%act = wx(j,inds(indPlot));
%         actI = max(1,min(64,32+round(10*real(act))));
%         %cols = squeeze(complexIm(act,0,.5));
%         subplot(5,5,j);scatter(real(Xor(inds(indPlot))),imag(Xor(inds(indPlot))),abs(act)*10,cols(actI,:),'filled');%hold all;
%         axis tight
%         drawnow;
%     end