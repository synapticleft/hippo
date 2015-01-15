function makeHists(fn,inds,window)
%% make 3 rows of 2-d histograms that show the joint distribution of two columns of data
%% each of 3 rows has all, left choices, right choices
%% each of 6 columns are for 6 ILD's.

if ~exist('figsOn','var')
    figsOn = 1;
end
if ~exist('window','var')
    window = [180 240];
end

file = dir('*.mat');
load(file(fn).name,'data');
file(fn).name

temp = [data{2:end,2}];

allData = nan*ones(numel(temp),window(2)-window(1),numel(inds));
trs = min(temp):max(temp);

for i = trs %1:max(temp)
    f = find(i == temp);
    for j = 1:numel(f)
        temp1 = reshape([data{f(j)+1,inds}],numel(data{f(j)+1,inds(1)}),[]);
        allData(f(j),:,:) = temp1(window(1)+1:window(2),:);
    end
end

%allData(isnan(allData)) = 0;
choice = [data{2:end,7}];
answer = [data{2:end,5}];
difficulty = [data{2:end,8}];
%difficulty = circshift(difficulty,[0 -1]);
f = choice == 3;% | abs(difficulty) > .5;
allData(f,:,:) = [];
choice(f) = [];
answer(f) = [];
temp(f) = [];
difficulty(f) = [];
binAnswer = double(answer > 0)*2 - 1;
choice = choice*2 - 3;

figure;
% correct the drift of eye from session to session, but not for other recordings
for j = 1:size(allData,3)
    if ismember(inds(j),10:14)
        for i = trs %1:max(temp)
            f = find(i == temp);
            temp1 = squeeze(allData(f,:,j));
            temp1 = temp1 - nanmean(temp1(:));
            temp1 = temp1 / nanstd(temp1(:));
            allData(f,:,j) = temp1;
        end
%    else
%        temp = squeeze(allData(:,:,j));
%        allData(:,:,j) = allData(:,:,j) - mean(temp(:));
%        allData(:,:,j) = allData(:,:,j) / std(temp(:));
    end
end

y = answer;%binAnswer;%choice;%

allData = permute(allData,[3 1 2]);
bins{1} = linspace(prctile(allData(1,:),1),prctile(allData(1,:),99),20);
bins{2} = linspace(prctile(allData(2,:),5),prctile(allData(2,:),95),20);
y = circshift(choice,[0 0]);%.*binAnswer
for i = 1:6
    f = difficulty == i;
    temp = allData(:,f,:);
    subplot(3,6,i);imagesc(bins{2},-bins{1},log(hist3(temp(:,:)',bins)));
    temp = allData(:,difficulty == i & y < 0,:);
    subplot(3,6,i+6);imagesc(bins{2},-bins{1},log(hist3(temp(:,:)',bins)));
    temp = allData(:,difficulty == i & y > 0,:);
    subplot(3,6,i+12);imagesc(bins{2},-bins{1},log(hist3(temp(:,:)',bins)));
end

%figure;imagesc(log(hist3(allData',bins)));