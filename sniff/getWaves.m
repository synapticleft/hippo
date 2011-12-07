function [Sniff lfp odors sniffTimes ] = getWaves(folder,bin,len,channels)%

f = dir([folder '/rsm*.mat']);
s = whos('-regexp','\d','-file',[folder '/' f(1).name]);
load([folder '/' f(1).name],'Sniff');%_DC
Sniff = decimate(Sniff(1:len),bin);%_DC
if nargout > 1
    if ~exist('channels','var')
        channels = 1:numel(s);
    end
    lfp = zeros(numel(channels),ceil(len/bin));%s(1).size(2)
    for i = 1:size(lfp,1)
        load([folder '/' f(1).name],s(i).name);
        temp = eval(s(i).name);
        lfp(i,:) = decimate(temp(1:len),bin);
    end
end
if nargout > 3
    f = dir([folder '/sniff*.mat']);
    load([folder '/' f(1).name],'sniff');
    temp = [sniff.t0] + sniff(1).t_zer(1);
    temp(temp > len) = [];
    sniffTimes = zeros(1,size(lfp,2));
    sniffTimes(round(temp/bin)) = 1;
end
if nargout > 2
    f = dir([folder,'/sess*.mat']);
    load([folder '/' f(1).name],'trial');
    odors = zeros(max([trial.odorValve]+1),size(lfp,2));
    for i = 1:numel(trial)
        times = trial(i).start + trial(i).odorTime;
        if times(2) > len
            break
        end
        times = round(times/bin);
        if diff(times) > 0
            try
                odors(trial(i).odorValve+1,times(1):times(2)) = 1;
            catch
                [trial(i).odorValve times']
            end
        end
    end    
end