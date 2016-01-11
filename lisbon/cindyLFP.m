function [lfpChunk,sniffChunk,data] = cindyLFP(lfp,sniff,meta,Fs,prePost)
%get snippets of data that correspond to periods where rat is sniffing, for
%LFP and sniff data. Also get a vector of origin, target (allo and ego),
%and RT.

valid = ~isnan(meta.behav_initmsd_on);%logical(meta.behav.init_valid);
start = meta.behav_initmsd_on/1000;%meta.behav.init_msd_delay(:,1);
%stop = meta.behav_goalpokein/1000;
%start = stop;

lfpChunk = zeros(size(lfp,1),sum(valid),diff(prePost)*Fs);
sniffChunk = zeros(sum(valid),diff(prePost)*Fs);

start = start(valid);

for i = 1:numel(start)
    lfpChunk(:,i,:) = lfp(:,round(Fs*start(i))+(prePost(1)*Fs+1:prePost(2)*Fs));
    sniffChunk(i,:) = sniff(round(Fs*start(i))+(prePost(1)*Fs+1:prePost(2)*Fs));
end
temp = meta.behav_goal;temp = circshift(temp,[0 1]);temp(1) = 0;
data = [meta.behav_init(valid)' meta.behav_goal(valid)' meta.trialstruct_goal(valid)' temp(valid)']; 
data1 = data;
data1(data == 2) = 3;
data1(data == 3) = 2;
data = data1;
data = [data  mod(data(:,1)-data(:,2),4)  meta.behav_initpokeout(valid)'-meta.behav_initmsd_on(valid)'];
lfpChunk = lfpChunk(:,2:end,:);sniffChunk = sniffChunk(2:end,:); data = data(2:end,:);
    % data = [meta.trialstruct.init(valid) meta.trialstruct.goal(valid) ...
%     mod(meta.trialstruct.init(valid)-meta.trialstruct.goal(valid),4) ...
%     meta.behav.init_pokeout(valid)-meta.behav.init_msd_delay(valid,1)];

% to do - classifier with sliding window; ripples; goal aligned; invalid
% trials