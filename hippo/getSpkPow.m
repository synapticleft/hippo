function allS = getSpkPow(file,shanks,v)
%% creates a power spectrum of spike waveforms for each neurons, then takes the hilbert transform in theta band

d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
nSpikes = inf;%20000;
params.Fs = 20000;
win = [32/1250 32/1250];
params.tapers = [1/win(1) win(1) 1];params.fpass = [500 3000];
%params1.Fs = 1250/32;params1.tapers = [15 29];
allS = zeros(numel(shanks)*8,size(v,1));counter = 0;
for i = shanks
    try
         spk = LoadSpk([d file '.spk.' num2str(i)],8,32,nSpikes);
    catch
        spk = LoadSpk([d file '.spk.' num2str(i)],7,32,nSpikes);
    end
     times = LoadResTimes([d file],i);
     times = times(1:min(nSpikes,numel(times)));
%     timeCourse = zeros(size(spk,1),32*16*size(v,1));
     for j = size(spk,2):-1:1
         timeCourse(:,j+times-1) = squeeze(spk(:,j,:));
     end
     if size(timeCourse,2) > size(v,1)*32*16
         timeCourse = timeCourse(:,1:size(v,1)*32*16);
     end
     [S,~,f] = mtspecgramc(timeCourse',win,params);
     clear timeCourse;
     S = squeeze(mean(S,2));
     Sf = morFilter(S',8,1250/32);
     if size(Sf,2) < size(v,1)
         Sf(:,size(Sf,2)+1:size(v,1)) = 0;
     end
     allS((counter+1):(counter+size(spk,1)),:) = Sf;
     clear Sf S;counter = counter + size(spk,1);
end