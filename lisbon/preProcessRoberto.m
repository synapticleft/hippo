function [allData, allOut, allOutShift,center_ILD] = preProcessRoberto(fn,inds,timePast,offSet,trs,whichDiff,startAlign,addClassic)
%% function used to organize all data into a matrix that can be used for further analysis.
%% inputs
%% fn: filename, used to select subject / stimulus condition
%% inds: which features to use for predictions
%% timePast: how many samples should be used to make prediction (used 60, or .75 s)
%% trs: sessions (default all)
%% whichDiff: which ILD's to use for analysis
%% startAlign: align to stimulus onset or offset.
%% addClassic: add subject's choice and mean ILD as regressors for prediction

window = [120 260];

file = dir('*.mat');
load(file(fn).name,'data','center_ILD');
file(fn).name
if ~exist('inds','var')
    inds = 10:15;%21;
end
if ~exist('whichDiff','var') || isempty(whichDiff)
    whichDiff = [-4 -2 -.5 .5 2 4];
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


f1 = find(ismember(inds,[10:14 15])); %%%ADDED HAND 
f = choice == 3 | ~ismember(difficulty,whichDiff);
if ~isempty(f1)
   f = f | squeeze(~(allData(:,1,f1(1))))';
   if fn == 1 || fn == 2
       f = f | [data{2:end,2}] < 4;
   end
end
allOut = [answer' double(answer' > 0)*2-1 choice'*2-3];
allOut = [allOut allOut(:,2).*allOut(:,3)];
allOut = [allOut  temp' difficulty'];
allOutShift = [circshift(allOut,[1 0]) circshift(allOut,[-1 0])];
%answer = circshift(answer,[0 1]);
%choice = circshift(choice,[0 1]);
allOut(f,:) = [];
allOutShift(f,:) = [];
allData(f,:,:) = [];
%choice(f) = [];
%answer(f) = [];
temp(f) = [];
%difficulty(f) = [];

%binAnswer = double(answer > 0)*2 - 1;
%choice = choice*2 - 3;
%correct = choice.*binAnswer;

% correct the drift of eye from session to session, but not for other recordings
for j = 1:numel(inds)%size(allData,3)
    if ismember(inds(j),10:14)
        for i = trs %1:max(temp)
            f = find(i == temp);
            temp1 = squeeze(allData(f,:,j));
            temp1 = temp1 - mean(temp1(:));
            temp1 = temp1 / std(temp1(:));
            allData(f,:,j) = temp1;
        end
    elseif ~ismember(inds(j),6)
       tempa = squeeze(allData(:,:,j));
       allData(:,:,j) = allData(:,:,j) - mean(tempa(:));
       allData(:,:,j) = allData(:,:,j) / std(tempa(:));
    end
end
if addClassic
    allData = circshift(allData,[0 0 2]);
end