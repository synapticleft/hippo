function [allData,center_ILD] = preProcessRoberto1(fn,inds,timePast,offSet,trs,whichDiff,startAlign,addClassic) %, allOut, allOutShift
%% function used to organize all data into a matrix that can be used for further analysis.
%% inputs
%% fn: filename, used to select subject / stimulus condition
%% inds: which features to use for predictions
%% timePast: how many samples should be used to make prediction (used 60, or .75 s)
%% trs: sessions (default all)
%% whichDiff: which ILD's to use for analysis
%% startAlign: align to stimulus onset or offset.
%% addClassic: add subject's choice and mean ILD as regressors for prediction

window = [120 300]; %was 260
maxNan = 5;
validSession = .3;

file = dir('*.mat');
load(file(fn).name,'data','center_ILD');
file(fn).name
if ~exist('inds','var')
    inds = 10:15;%21;
end
if ~exist('whichDiff','var') || isempty(whichDiff)
    whichDiff = [-6 -4 -2 -.5 .5 2 4 6];
end
if ~exist('timePast','var')
    timePast = 0;
end
if ~exist('offSet','var')
    offSet = 0;
end
if ~exist('startAlign','var')
    startAlign = 1;
end
if ~exist('addClassic','var')
    addClassic = 0;
end

temp = [data{2:end,2}];
if ~exist('trs','var') || isempty(trs)
    trs = min(temp):max(temp);
end

allData = nan*ones(numel(temp),window(2)-window(1)+timePast,numel(inds));


for i = trs 
    f = find(i == temp);
    for j = 1:numel(f)
        temp1 = reshape([data{f(j)+1,inds}],numel(data{f(j)+1,inds(1)}),[]);
        if startAlign
            allData(f(j),:,1:end-1) = temp1((window(1)-timePast+1:window(2))-offSet,1:end-1);
            allData(f(j),:,end) = temp1((window(1)-timePast+1:window(2)),end);
        else
            allData(f(j),:,1:end-1) = temp1(end-(window(2)+timePast-1:-1:window(1))-offSet,1:end-1);
            allData(f(j),:,end) = temp1(end-(window(2)+timePast-1:-1:window(1)),end);
        end
    end
end
if addClassic
    %allData(:,:,end+1) = 0;
    %for i = trs
    %    f = find(i == temp);
    %    t = [data{f+1,9}]';
    %    t = t-nanmean(t);t = t/nanstd(t);
    %    allData(f,:,end) = repmat(t,[1 size(allData,2) ]);
    %end
    allData(:,:,end+1) = repmat(([data{2:end,7}]'-1.5),[1 size(allData,2) 1]);
    allData(:,:,end+1) = repmat(([data{2:end,5}]' > center_ILD),[1 size(allData,2) 1]);
    %allData(:,:,end+1) = allData(:,:,end-1).*allData(:,:,end);%repmat([data{2:end,9}]',[1 size(allData,2) 1]);
end
f = inds == 6;
if startAlign
    allData(:,nanmean(abs(diff(allData(:,:,f),[],2))) == 0,f) = nan;
end
%allData(:,:,f) = bsxfun(@minus,allData(:,:,f),nanmean(allData(:,:,f),2));
%allData(:,:,f) = allData(:,:,f) + eps;
allData(isnan(allData)) = 0;%(:,1,f)),:,f) = 0;
choice = [data{2:end,7}];
answer = [data{2:end,5}];
difficulty = [data{2:end,4}];


f1 = find(ismember(inds,[10:14])); %%%15 ADDED HAND 
f = choice == 3;
%figure;
if ~isempty(f1)
    for i = 2:size(data,1)
        tempa = bwlabel(isnan(data{i,24}));
        ts = [find(~isnan(data{i,6}),1) find(~isnan(data{i,6}),1,'last')];
        %ts = window(1) + [find(allData(i-1,timePast + (1:60),end) == 0,1,'last')+1 60+find(allData(i-1,timePast + (61:end),end) == 0,1)+1];
        %if numel(ts) == 1
        %    ts(2) = window(1) + size(allData,2);
        %end
        %plot(ts(1):ts(2),data{i,6}(ts(1):ts(2)));input('');
        ts = data{i,22}(ts);
        temp1 = data{i,28} > ts(1) & data{i,28} < ts(2);
        temp1 = tempa.*temp1';
        for j = min(temp1(temp1 > 0)):max(temp1)
            f(i-1) = f(i-1) | sum(tempa == j) > maxNan;
        end
%         if f(i-1) && choice(i-1) ~= 3
%             hold all;plot(data{i,23}(data{i,28} > ts(1) & data{i,28} < ts(2)));
%         end
    end
end
running = zeros(numel(f),1+numel(f1));
for i = find(~f)
    running(i,1) = sum(~isnan((data{i+1,23})));
    for j = 1:numel(f1)
        running(i,j+1) = nanmean(data{i+1,f1(j)+22});
        running(i,j+4) = nanstd(data{i+1,f1(j)+22});
    end
end
running(running == 0) = nan;
%figure;plot(running(:,end-1:end));
%return
for i = 1:max(temp)
    if sum(~f(temp == i))/sum(temp == i) < validSession
        %sum(~f(temp == i))/sum(temp == i)
        f(temp == i) = 1;
    end
end
 f = f | ~ismember(difficulty,whichDiff);
% hold all;plot(f,'r')
%[sum(f & choice ~= 3) sum(f) sum(choice == 3)]/numel(f)

% allOut = [answer' double(answer' > 0)*2-1 choice'*2-3];
% allOut = [allOut allOut(:,2).*allOut(:,3)];
% allOut = [allOut  temp' difficulty'];
% allOutShift = [circshift(allOut,[1 0]) circshift(allOut,[-1 0])];
% %answer = circshift(answer,[0 1]);
% %choice = circshift(choice,[0 1]);
% allOut(f,:) = [];
% allOutShift(f,:) = [];
%figure;plot(f);hold all;plot(choice == 3,'r');
%choice(f) = [];
%answer(f) = [];
%difficulty(f) = [];
%binAnswer = double(answer > 0)*2 - 1;
%choice = choice*2 - 3;
%correct = choice.*binAnswer;
%figure;
% correct the drift of eye from session to session, but not for other recordings
for j = 1:numel(inds)%size(allData,3)
    if ismember(inds(j),10:14)
        for i = trs %1:max(temp)
            f1 = find(i == temp);
            if 1
               % ff = find(f1);
                if sum(i == temp & ~f)%~isempty(ff)
                    temp1 = interp1(find(i == temp & ~f),running(i == temp & ~f,1+inds(j)-9),f1,'pchip');
                    %temp1 = filtfilt(gausswin(10),sum(gausswin(10)),temp1);
                    temp1 = smooth(temp1,19,'rloess');
                    allData(f1,:,j) = bsxfun(@minus,allData(f1,:,j),temp1);
                    temp1 = squeeze(allData(f1,:,j));
                     temp1 = temp1 / std(temp1(:));
                    allData(f1,:,j) = temp1;
%                    plot(f1,temp1);hold all;plot(find(i == temp & ~f),running(i == temp & ~f,1+inds(j)-9));hold off; pause(.5);
                end
                %%SMOOTH TO REMOVE IDIOSYNCRACIES
            else
                temp1 = squeeze(allData(f1,:,j));
                temp1 = temp1 - mean(temp1(:));
                temp1 = temp1 / std(temp1(:));
                allData(f1,:,j) = temp1;
            end
        end
    elseif ~ismember(inds(j),6)
       tempa = squeeze(allData(:,:,j));
       allData(:,:,j) = allData(:,:,j) - mean(tempa(:));
       allData(:,:,j) = allData(:,:,j) / std(tempa(:));
    end
end
allData(f,:,:) = [];
% figure;subplot(211);plot(squeeze(mean(allData(:,:,1),2)));
% hold all;plot(squeeze(mean(allData(:,:,2),2)));
% subplot(212);scatter(squeeze(allData(:,40,1)),squeeze(allData(:,100,1)));
% hold all;scatter(squeeze(allData(:,40,2)),squeeze(allData(:,100,2)));
% temp(f) = [];
if addClassic
    allData = circshift(allData,[0 0 2]);
end