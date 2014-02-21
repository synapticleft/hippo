function [spikeMat cellInfo] = hipSpikes(file,bin,subSet,d)%spk %a b c badUnits
%% extract all spikes into matrix downsampled to match lfp and position files.
%bin usually 32/1.25
%%new way for Kenji's data

%load('/media/Expansion Drive/redwood/KenjiData.mat');
load('/media/work/hippocampus/KenjiData.mat');
whichDay = strcmp(file,Beh(:,4));
dayID = Beh(whichDay,2);
dayCells = PyrIntMap.Map(:,1) == find(strcmp(dayID,PyrIntMap.fileBase)) & Region  == 1;%> 3;%
someShanks = unique(PyrIntMap.Map(dayCells,3));
%%old way
if ~exist('d','var') || isempty(d)
    %d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
    d = ['/media/Kenji_data/' Beh{whichDay,3} '/' Beh{whichDay,1} '/' file '/'];
end
if exist('someShanks','var') && ~isempty(someShanks)
    [a,b,c,d] = LoadCluRes([d file],someShanks);
else
     [a,b,c,d] = LoadCluRes([d file]);
end
a = ceil(a/d.SampleRate*1000/bin);
if ~exist('subSet','var')
    subSet = [0 max(a)];
else
%    subSet = subSet*1000/bin;
end

b(a < subSet(1)) = [];a(a < subSet(1)) = [];
if numel(subSet) == 2
    b(a > subSet(2)) = [];
    a(a > subSet(2)) = [];
end
a = a - subSet(1);
spikeMat = zeros(sum(dayCells),max(a));
cellInfo.type = -20*ones(sum(dayCells),1);
cellInfo.shank = cellInfo.type;
cellInfo.ID = cellInfo.type;
cellInfo.popSize = sum(dayCells);
for i = 1:max(b)
    temp = hist(a(b==i),1:(max(a)+1));
    f = find(c(i,2) == PyrIntMap.Map(dayCells,3) & c(i,3) == PyrIntMap.Map(dayCells,4));
    if numel(f)
    assert(numel(f) == 1);
    spikeMat(f,:) = temp(1:(end-1));%a(b == i)) = 1;
    if c(i,3) == 0 || c(i,3) == 1
        cellInfo.type(f) = c(i,3) - 2;
    else
        cellInd = dayCells & c(i,2) == PyrIntMap.Map(:,3) & c(i,3) == PyrIntMap.Map(:,4);
        cellInfo.ID(f) = PyrIntMap.Map(cellInd,2);
        cellInfo.shank(f) = c(i,2);
        assert(sum(cellInd) == 1);
        if ~Clean(cellInd)
            cellInfo.type(f) = 0;
        else
            for j = 1:3
                if iCell{j}(cellInd)
                    cellInfo.type(f) = j;
                end
            end
        end
    end
    end
end
%badUnits = [ spikeMat(c(:,3) == 0,:); spikeMat(c(:,3) == 1,:)];
%spikeMat(c(:,3) ==0 | c(:,3) == 1,:) = [];
