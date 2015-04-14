function [allAct allActC allStim allStimC myShift] = visICAAct(act,trials,stim,actInd)%pos,

allAct = [];
trials(trials == 0) = [];
trials = [0 trials];
for i = 1:numel(trials)-1
allAct(:,1:(trials(i+1)-trials(i)),i) = act(:,trials(i)+1:trials(i+1));
allStim(1:(trials(i+1)-trials(i)),i) = stim(trials(i)+1:trials(i+1));
end


figure;showGrid(permute(allAct,[1 3 2]));
if ~exist('actInd','var')
    actInd = input('Which component is reference? ');
end

refAct = (squeeze(allAct(actInd,:,:)));

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
for i = 1:size(refAct,2)
    %myShift(i) = find(abs(filtfilt(gausswin(5),sum(gausswin(5)),double(abs(filtfilt(gausswin(5),sum(gausswin(5)),refAct(:,i))) > .5)) - 1) < .01,1);
    temp(i,:) = abs(filtfilt(gausswin(5),sum(gausswin(5)),double(abs(filtfilt(gausswin(5),sum(gausswin(5)),refAct(:,i))) > .5)) - 1) < .01;
    temp(i,1) = 0;myShift(i) = find(temp(i,:),1);
%    b = bwlabel(f);b = hist(b,0:max(b));b = b(2:end);
    %refActC(:,i) = circshift(refAct(:,i),[-myShift(i) 0]);
end

myShift = myShift - round(mean(myShift));
allStimC = allStim;
for i = 1:size(refAct,2)
    allActC(:,:,i) = circshift(allAct(:,:,i),[0 -myShift(i) 0]);
    allStimC(:,i) = circshift(allStim(:,i),[-myShift(i) 0]);
end

%    [~,m] = max(conv([zeros(chunk,1); refActC(:,i)],flipud(refActm)));
%    
% end

%figure;showGrid(permute(allActC,[1 3 2]));%subplot(211);imagesc(abs(squeeze(allAct(actInd,:,:))'));subplot(212);imagesc(abs(squeeze(allActC(actInd,:,:))'));

%figure;imagesc(temp);
[~,s] = sort(mean(allStim(1:20,:)));
allStim = allStim(:,s);
allStimC = allStimC(:,s);
figure;showGrid(permute(allAct(:,:,s),[1 3 2]));
figure;showGrid(permute(allActC(:,:,s),[1 3 2]));
allActC = allActC(:,:,s);
allAct = allAct(:,:,s);
% allHist = [];
% g = meshgrid(1:size(allAct,2),1:size(allAct,3))';
% xs{1} = -3:.1:3;xs{2} = 1:size(allAct,2);
% for i = 1:size(allAct,1)
%     t = real(allActC(i,:)) ~= 0;
%     allHist(i,:,:) = hist3([real(allActC(i,t))' g(t)'],xs);
% end
% figure;showGrid(allHist);
