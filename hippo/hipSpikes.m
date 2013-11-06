function [spikeMat cellInfo] = hipSpikes(file,bin,d,subSet)%spk %a b c badUnits
%% extract all spikes into matrix downsampled to match lfp and position files.
%bin usually 32/1.25
%%new way for Kenji's data

%load('/media/Expansion Drive/redwood/KenjiData.mat');
load('/media/work/hippocampus/KenjiData.mat');
whichDay = strcmp(file,Beh(:,4));
dayID = Beh(whichDay,2);
dayCells = PyrIntMap.Map(:,1) == find(strcmp(dayID,PyrIntMap.fileBase)) & Region  == 1;%> 3;%
%d = ['/media/Kenji_data/' Beh{whichDay,3} '/' Beh{whichDay,1} '/' file '/'];
someShanks = unique(PyrIntMap.Map(dayCells,3));
%%old way
if ~exist('d','var') || isempty(d)
    d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
end
if exist('someShanks','var') && ~isempty(someShanks)
    [a,b,c,d] = LoadCluRes([d file],someShanks);
else
     [a,b,c,d] = LoadCluRes([d file]);
end
a = ceil(a/d.SampleRate*1000/bin);
if ~exist('subSet','var')
    subSet = max(a);
else
    subSet = subSet*1000/bin;
end
mb = max(b);
b(a > subSet) = [];a(a > subSet) = [];
spikeMat = zeros(max(b),max(a));
%cellInfo = -20*ones(1,max(b));
%shankInfo = c(:,2);
cellInfo.type = -20*ones(size(c,1),1);
cellInfo.shank = c(:,2);
cellInfo.ID = cellInfo.type;
cellInfo.popSize = sum(dayCells);
for i = 1:mb%max(b)
    temp = hist(a(b==i),1:(max(a)+1));
    spikeMat(i,:) = temp(1:(end-1));%a(b == i)) = 1;
    if c(i,3) == 0 || c(i,3) == 1
        cellInfo.type(i) = c(i,3) - 2;
    else
        cellInd = dayCells & c(i,2) == PyrIntMap.Map(:,3) & c(i,3) == PyrIntMap.Map(:,4);
        cellInfo.ID(i) = PyrIntMap.Map(cellInd,2);
        assert(sum(cellInd) == 1);
        if ~Clean(cellInd)
            cellInfo.type(i) = 0;
        else
            for j = 1:3
                if iCell{j}(cellInd)
                    cellInfo.type(i) = j;
                end
            end
        end
    end
end
%badUnits = [ spikeMat(c(:,3) == 0,:); spikeMat(c(:,3) == 1,:)];
%spikeMat(c(:,3) ==0 | c(:,3) == 1,:) = [];
