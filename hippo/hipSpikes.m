function [spikeMat a b c badUnits] = hipSpikes(file,bin,someShanks,subSet)
%% extract all spikes into matrix downsampled to match lfp and position files.
%bin usually 32/1.25
d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
d = '/media/Kenji_data/unit3/ec014.29/ec014.468/';
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
for i = 1:max(b)
    temp = hist(a(b==i),1:(max(a)+1));
    spikeMat(i,:) = temp(1:(end-1));%a(b == i)) = 1;
end
size(spikeMat)
badUnits = [ spikeMat(c(:,3) == 0,:); spikeMat(c(:,3) == 1,:)];
spikeMat(c(:,3) ==0 | c(:,3) == 1,:) = [];
