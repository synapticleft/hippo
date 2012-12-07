function spk = hipSpikes(file,bin,subSet)%[spikeMat cellInfo shankInfo]%a b c badUnits
%% extract all spikes into matrix downsampled to match lfp and position files.
%bin usually 32/1.25
%%new way for Kenji's data
load('/media/work/hippocampus/KenjiData.mat');
whichDay = strcmp(file,Beh(:,4));
dayID = Beh(whichDay,2);
dayCells = PyrIntMap.Map(:,1) == find(strcmp(dayID,PyrIntMap.fileBase)) & Region == 1;
d = ['/media/Kenji_data/' Beh{whichDay,3} '/' Beh{whichDay,1} '/' file '/'];
someShanks = unique(PyrIntMap.Map(dayCells,3));
%%old way
%d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
for i = 1:8
spk{i} = LoadSpk([d file '.spk.' num2str(i)],8);
end
return
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
b(a > subSet) = [];a(a > subSet) = [];
spikeMat = zeros(max(b),max(a));
cellInfo = -20*ones(1,max(b));
shankInfo = c(:,2);
for i = 1:max(b)
    temp = hist(a(b==i),1:(max(a)+1));
    spikeMat(i,:) = temp(1:(end-1));%a(b == i)) = 1;
    if c(i,3) == 0 || c(i,3) == 1
        cellInfo(i) = c(i,3) - 2;
    else
        cellInd = dayCells & c(i,2) == PyrIntMap.Map(:,3) & c(i,3) == PyrIntMap.Map(:,4);
        assert(sum(cellInd) == 1);
        if ~Clean(cellInd)
            cellInfo(i) = 0;
        else
            for j = 1:3
                if iCell{j}(cellInd)
                    cellInfo(i) = j;
                end
            end
        end
    end
end
%badUnits = [ spikeMat(c(:,3) == 0,:); spikeMat(c(:,3) == 1,:)];
%spikeMat(c(:,3) ==0 | c(:,3) == 1,:) = [];
