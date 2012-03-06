function [spikeMat c] = hipSpikes(file,bin,subSet)
%% extract all spikes into matrix downsampled to match lfp and position files.
%bin usually 32/1.25
d = '/media/Expansion Drive/KenjiMizuseki/';%['/media/work/hippocampus/' file '/'];%
[a,b,c,d] = LoadCluRes([d file]);
a = ceil(a/d.SampleRate*1000/bin);
if ~exist('subSet','var')
    subSet = max(a);
else
    subSet = subSet*1000/bin;
end
subSet
b(a > subSet) = [];a(a > subSet) = [];
spikeMat = zeros(max(b),max(a));
for i = 1:max(b)
    temp = hist(a(b==i),1:(max(a)+1));
    spikeMat(i,:) = temp(1:(end-1));%a(b == i)) = 1;
end

