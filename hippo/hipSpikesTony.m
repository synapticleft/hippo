function [spikeMat cellInfo] = hipSpikesTony(file,dec,subSet)%a b c badUnits
%% extract all spikes into matrix downsampled to match lfp and position files.
%% for cells not in Kenji's metadata
%bin usually 32/1.25
%%new way for Kenji's data
%load('/media/work/hippocampus/KenjiData.mat');
%whichDay = strcmp(file,Beh(:,4));
%dayID = Beh(whichDay,2);
%dayCells = PyrIntMap.Map(:,1) == find(strcmp(dayID,PyrIntMap.fileBase));% & Region == 1;
%d = ['/media/Kenji_data/' Beh{whichDay,3} '/' Beh{whichDay,1} '/' file '/'];
%someShanks = unique(PyrIntMap.Map(dayCells,3));
%%old way
%d = ['/media/work/hippocampus/AB3-58/'];% file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
%if exist('someShanks','var') && ~isempty(someShanks)
%    [a,b,c,d] = LoadCluRes([d file],someShanks);
%else
    [a,b,c,d] = LoadCluRes([file]);
%end
a = ceil(a/d.SampleRate*1250);
if ~exist('subSet','var')
    subSet = max(a);
else
    %subSet = subSet*1000/bin;
end
b(a <= subSet(1)*32) = [];a(a <= subSet(1)*32) = [];
b(a > subSet(2)*32) = [];a(a > subSet(2)*32) = [];
a = a-subSet(1)*32;
cellInfo = -20*ones(1,max(b));
for i = max(b):-1:1
    tic;
%    temp = hist(a(b==i),1:(max(a)+1));
%    spikeMat(i,:) = temp(1:(end-1));
    temp = hist(a(b == i),1:(diff(subSet)*32+1));
    spikeMat(i,:) = decimate(temp(1:(end-1)),dec);
    if c(i,3) == 0 || c(i,3) == 1
        cellInfo(i) = c(i,3) - 2;
    else
        cellInfo(i) = 1;
    end
    toc
end
