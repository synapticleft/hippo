function [coeff,cc] = evolvePredictTraj(fn,lambda,inds,timePast,offSet,trialHist,figsOn) %sessNorm,
% fit LDA and regularized linear regression of start- and end- aligned
% trials, look at performance, and weights.

% choices: 1) re-normalize measures in each session to correct for drifts
% (sessNorm) NOW I DO THIS FOR EYE DATA BUT NOT OTHERS
% 2) regularization of fit (lambda)
% 3) fit choice, previous choice, or right answer
% 4) number of time steps in the past
% 5) which measurements to include
% 6) separate each difficulty level or combine them all

if ~exist('timePast','var')
    timePast = 30;
end
if ~exist('figsOn','var')
    figsOn = 1;
end
if ~exist('trialHist','var')
    trialHist = 0;
end
window = [120 300];

file = dir('*.mat');
load(file(fn).name,'data');
file(fn).name
if ~exist('inds','var')
    inds = 10:15;%21;
end
temp = [data{2:end,2}];

allData = nan*ones(numel(temp),window(2)-window(1)+timePast,numel(inds));

trs = min(temp):max(temp);

for i = trs %1:max(temp)
    f = find(i == temp);
    for j = 1:numel(f)
        temp1 = reshape([data{f(j)+1,inds}],numel(data{f(j)+1,inds(1)}),[]);
        allData(f(j),:,1:end-1) = temp1(window(1)-timePast+1:window(2),1:end-1);
        allData(f(j),:,end) = temp1(offSet+(window(1)-timePast+1:window(2)),end);
    end
end
f = inds == 6;
allData(:,mean(abs(diff(allData(:,:,1),[],2))) == 0,f) = nan;
allData(isnan(allData)) = 0;%(:,1,f)),:,f) = 0;
choice = [data{2:end,7}];
answer = [data{2:end,5}];
difficulty = [data{2:end,4}];


f1 = find(ismember(inds,10:14));
f = choice == 3 | difficulty ~= 4;%abs(difficulty) > .5;%abs(difficulty) > .5;
if ~isempty(f1)
   f = f | squeeze(~(allData(:,1,f1(1))))';
end
answer = circshift(answer,[0 1]);
choice = circshift(choice,[0 1]);
allData(f,:,:) = [];
choice(f) = [];
answer(f) = [];
temp(f) = [];
difficulty(f) = [];
% f = abs(difficulty) > .5;
% allData(f,:,:) = [];
% choice(f) = [];
% answer(f) = [];
% temp(f) = [];
binAnswer = double(answer > 0)*2 - 1;
choice = choice*2 - 3;
correct = choice.*binAnswer;


% correct the drift of eye from session to session, but not for other recordings
for j = 1:size(allData,3)
    if ismember(inds(j),10:14)
        for i = trs %1:max(temp)
            f = find(i == temp);
            temp1 = squeeze(allData(f,:,j));
            temp1 = temp1 - mean(temp1(:));
            temp1 = temp1 / std(temp1(:));
            allData(f,:,j) = temp1;
        end
%    else
%        temp = squeeze(allData(:,:,j));
%        allData(:,:,j) = allData(:,:,j) - mean(temp(:));
%        allData(:,:,j) = allData(:,:,j) / std(temp(:));
    end
end

%y = answer;%binAnswer;%choice;%
%y = circshift(y,[0 1]);
%figure;hold all;
%allData(:,:,end) = cumsum(allData(:,:,end),2);
for i = timePast+1:size(allData,2)
    X = allData(:,i-timePast+1:i,1:end-1);
    y = zscore(squeeze(allData(:,i,end)))';
    X = X(:,:);
    if trialHist
        X = [X circshift([answer]',[0 0])];% circshift([answer; binAnswer; choice;correct]',[-1 0])];
    end
    X = zscore(X);
    coeff(i-timePast,:) = (X'*y')'/(X'*X + lambda*size(X,1)*eye(size(X,2)));%y/X';
    yHat = coeff(i-timePast,:)*X';
    mse(i-timePast) = mean(sqrt((y-yHat).^2));%sum(yHat'+1 ~=info(inds,class))/sum(inds);
    temp = corrcoef(y,yHat);
    cc(i-timePast) = temp(1,2);
    if 0 % figsOn
    if (i-timePast == 90) || (i-timePast == 120)
        scatter(y,yHat,'filled');
    end
    end
end
if figsOn
%coeff = reshape(coeff,size(coeff,1),timePast,[]);
figure;imagesc(coeff(:,:)');
%figure;plot(coeff(:,:));
figure;plot(max(mse)-mse);hold all;plot(cc);ylim([0 1]);
[~,m] = min(mse);
%figure;plot(squeeze(coeff(90,:,:)));
end