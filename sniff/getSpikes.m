function [sp ] = getSpikes(folder,bin,len,cells)

f = dir([folder '/spikes*.mat']);
load([folder '/' f(1).name],'unit');
if ~exist('cells','var')
    cells = 1:numel(unit);
end
if ~exist('len','var')
    len = 0;
    for i = 1:numel(cells)
        len = max(len,ceil(max(unit(cells(i)).spikeTimes)));
    end
end

sp = zeros(numel(cells),ceil(len/bin));

for i = 1:size(sp,1)
    times = unit(cells(i)).spikeTimes;
    times = ceil(times/bin);
    times(times > size(sp,2)) = [];
    sp(i,:) = hist(times,1:size(sp,2));
end