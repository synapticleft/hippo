function [allAct allActC allStim allStimC myShift] = visICAAct(act,trials,stim,refInd,stimInd)%pos,
%% do end-align

startAlign = 1;
stimRefVal = stim(trials(1:end-1) + 14);
allAct = [];
%trials(trials == 0) = [];
%trials = [0 trials];
for i = 1:numel(trials)-1
allAct(:,1:(trials(i+1)-trials(i)),i) = act(:,trials(i)+1:trials(i+1));
allStim(1:(trials(i+1)-trials(i)),i) = stim(trials(i)+1:trials(i+1));
end
    if ~startAlign
        sz = diff(trials);
        for i = 1:size(allAct,3)
            allAct(:,:,i) = circshift(allAct(:,:,i),[0 size(allAct,2)-sz(i) 0]);
            allStim(:,i) = circshift(allStim(:,i),[size(allAct,2)-sz(i) 0]);
        end
        allAct = allAct(:,end:-1:1,:);
        allStim = allStim(end:-1:1,:);
    end

%
if ~exist('refInd','var')  || isempty(refInd)
    figure;showGrid(permute(allAct,[1 3 2]));
    refInd = input('Which component is reference? ');
end
%return
refAct = (squeeze(allAct(refInd,:,:)));

% chunk = 40;
% refActm = mean(refAct,2);
% refActm = refActm - mean(refActm);
% refActm = refActm(1:chunk);
% %refActm = [-ones(1,chunk/2) ones(1,chunk/2)];
% for i = 1:size(refAct,2)
%     [~,m] = max(conv([zeros(chunk,1); refAct(:,i)],flipud(refActm)));
%     refActC(:,i) = circshift(refAct(:,i),[-m-chunk-20 0]);
% end
% refActm = mean(refActC,2);
% refActm = refActm(1:chunk);
allActC = allAct;
temp = [];
if startAlign
for i = 1:size(refAct,2)
    %myShift(i) = find(abs(filtfilt(gausswin(5),sum(gausswin(5)),double(abs(filtfilt(gausswin(5),sum(gausswin(5)),refAct(:,i))) > .5)) - 1) < .01,1);
    temp(i,:) = abs(filtfilt(gausswin(5),sum(gausswin(5)),double(abs(filtfilt(gausswin(5),sum(gausswin(5)),refAct(:,i))) > .5)) - 1) < .01;
    temp(i,1) = 0;
    myShift(i) = find(temp(i,:),1);
%    b = bwlabel(f);b = hist(b,0:max(b));b = b(2:end);
    %refActC(:,i) = circshift(refAct(:,i),[-myShift(i) 0]);
end
else 
    myShift = zeros(1,size(refAct,2));
end

myShift = myShift - round(median(myShift));
allStimC = allStim;
for i = 1:size(refAct,2)
    allActC(:,:,i) = circshift(allAct(:,:,i),[0 -myShift(i) 0]);
    allStimC(:,i) = circshift(allStim(:,i),[-myShift(i) 0]);
end

[~,s] = sort(mean(allStim(1:20,:)));
m1 = mdscale(squareform(pdist(squeeze(real(allActC(stimInd,21:60,:)))')),1);
[~,ms] = sort(m1);
figure;for i = 1:size(act,1)
subtightplot(4,size(act,1),i);imagesc(real(squeeze(allAct(i,1:100,:)))',[-2 2]);colormap jet;axis off;
subtightplot(4,size(act,1),i+size(act,1));imagesc(real(squeeze(allActC(i,1:100,:)))',[-2 2]);colormap jet;axis off;
subtightplot(4,size(act,1),i+size(act,1)*2);imagesc(real(squeeze(allActC(i,1:100,ms)))',[-2 2]);colormap jet;axis off;
subtightplot(4,size(act,1),i+size(act,1)*3);imagesc(real(squeeze(allActC(i,1:100,s)))',[-2 2]);colormap jet;axis off;
end

%figure;showGrid(permute(allActC,[1 3 2]));

allStim = allStim(:,s);
allStimC = allStimC(:,s);
allActC = allActC(:,:,s);
allAct = allAct(:,:,s);
stimRefVal = stimRefVal(s);

%figure;showGrid(permute(allAct,[1 3 2]));
%figure;showGrid(permute(allActC(:,1:100,:),[1 3 2]));

inds(1,:) = abs(stimRefVal) < 1;
inds(2,:) = abs(stimRefVal) < 3 & ~inds(1,:);
inds(3,:) = ~inds(2,:)  & ~inds(1,:);
sgn = -1;allStim = allStim*sgn;stimRefVal = stimRefVal*sgn;allStimC = allStimC*sgn;
if ~exist('stimInd','var') || isempty(stimInd)
stimInd = input('Which component is stim-sensitive? ');
end
stimAct = squeeze(allAct(stimInd,:,:));
figure; for i = 1:3
    subplot(1,3,i);imagesc(corr(stimAct(1:100,inds(i,:))',allStim(1:100,inds(i,:))'),[-0.05 .3]);axis square;
end
stimActC = squeeze(allActC(stimInd,:,:));
figure; for i = 1:3
    subplot(1,3,i);imagesc(corr(stimActC(1:100,inds(i,:))',allStim(1:100,inds(i,:))'),[-0.05 .3]);axis square;
end

% figure;subplot(221);imagesc(corr(stimAct(1:100,:)',allStim(1:100,:)'),[-.1 .3]);axis square
% xc = (stimAct(1:100,:)*stimAct(1:100,:)');
% xca = (allStim(1:100,:)*allStim(1:100,:)');
% lambda = 2;lambda1 = 10;
% subplot(222);imagesc((xc + lambda1*diag(diag(xc)))\stimAct(1:100,:)*1000*allStim(1:100,:)'/(xca + lambda1*diag(diag(xca))));
% subplot(223);imagesc((xc + lambda*diag(diag(xc)))\stimAct(1:100,:)*allStim(1:100,:)');
% subplot(224);imagesc(stimAct(1:100,:)*allStim(1:100,:)'/(xca + lambda*diag(diag(xca))));

% figure; for i = 1:3
%     xc = (stimAct(1:100,inds(i,:))*stimAct(1:100,inds(i,:))');
%     xca = (allStim(14:100,inds(i,:))*allStim(14:100,inds(i,:))');
%     subplot(1,3,i);imagesc((xc + lambda1*diag(diag(xc)))\stimAct(1:100,inds(i,:))...
%         *1000*allStim(14:100,inds(i,:))'/(xca + lambda1*diag(diag(xca))),[-0.05 .3]);axis square;
% end
% return
stims = [0 1 3 5 inf];%.5 2 4 inf];
figure;hist(stimRefVal,-15:.1:15);
%stims = [-.1 1 3 5 7 inf];
% figure;
% for i = 2:numel(stims)
%     for j = 1:2
%         inds = abs(allStim(14,:)) < stims(i) & abs(allStim(14,:)) > stims(i-1);
%         %inds = abs(allStim(14,:) - stims(i)*(-1)^j) < .1;
%         %sum(inds)%numel(stims)*(j-1)+i);
%         subplot(2,numel(stims)-1,i-1);imagesc(corr(stimAct(1:100,inds & allStim(14,:) > 0)',allStim(1:100,inds & allStim(14,:) > 0)'),[-.1 .3]);axis square;
%         title(sum(inds & allStim(14,:) > 0));
%         subplot(2,numel(stims)-1,i+numel(stims)-2);imagesc(corr(stimAct(1:100,inds & allStim(14,:) <= 0)',allStim(1:100,inds & allStim(14,:) <= 0)'),[-0.1 .3]);axis square;
%         title(sum(inds & allStim(14,:) <= 0));
%     end
% end
figure;
for i = 2:numel(stims)
        inds = abs(stimRefVal) < stims(i) & abs(stimRefVal) > stims(i-1);
        ccs(:,:,i-1) = (corr(stimAct(1:100,inds & stimRefVal > 0)',allStim(1:100,inds & stimRefVal > 0)')...
            +corr(stimAct(1:100,inds & stimRefVal <= 0)',allStim(1:100,inds & stimRefVal <= 0)'))/2;
        subplot(1,numel(stims)-1,i-1);imagesc(squeeze(ccs(:,:,i-1)),[-0.1 .3]);axis square;
        title(sum(inds));colormap jet
end
ccs(:,:,1) = mean(ccs(:,:,[1 2]),3);
ccs(:,:,2) = mean(ccs(:,:,[3 4]),3);ccs(:,:,[3 4]) = [];

for i = 1:100
ccs1(:,i,:) = circshift(ccs(:,i,:),[20-i 0 0]);
end
figure;for i = 1:2
    subplot(2,2,i);imagesc(ccs(:,:,i),[-.1 .3]);axis square;colormap jet
subplot(2,2,i+2);imagesc(ccs1(:,:,i),[-.1 .3]);axis square;colormap jet
end
figure;imagesc(corr(stimAct(1:100,:)',bsxfun(@minus,allStim(15:100,:),stimRefVal)'));

%doesnt work because groups multiple ILDs together in weird ways
% numdivs = 4;meanAct = mean(stimAct(1:30,:));
% figure;
% for i = 1:numdivs
%    subplot(2,2,i);inds = meanAct < prctile(meanAct,100*i/numdivs) & meanAct > prctile(meanAct,100*(i-1)/numdivs);
%    imagesc(corr(stimAct(1:100,inds)',allStim(1:100,inds)'-mean(mean(allStim(14:80,inds)))),[-.1 .3]);
% end
% allHist = [];
% g = meshgrid(1:size(allAct,2),1:size(allAct,3))';
% xs{1} = -3:.1:3;xs{2} = 1:size(allAct,2);
% for i = 1:size(allAct,1)
%     t = real(allActC(i,:)) ~= 0;
%     allHist(i,:,:) = hist3([real(allActC(i,t))' g(t)'],xs);
% end
% figure;showGrid(allHist);
