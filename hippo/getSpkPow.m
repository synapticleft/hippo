function Xp = getSpkPow(file,shanks,v)

d = ['/media/work/hippocampus/' file '/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
nSpikes = inf;%20000;
params.Fs = 20000;
win = [32/1250 32/1250];
params.tapers = [1/win(1) win(1) 1];
params1.Fs = 1250/32;params1.tapers = [15 29];
for i = shanks
    try
         spk = LoadSpk([d file '.spk.' num2str(i)],8,32,nSpikes);
    catch
        spk = LoadSpk([d file '.spk.' num2str(i)],7,32,nSpikes);
    end
     times = LoadResTimes([d file],i);
     times = times(1:min(nSpikes,numel(times)));
     timeCourse = zeros(size(spk,1),32*16*size(v,1));
     for j = size(spk,2):-1:1
         timeCourse(:,j+times-1) = squeeze(spk(:,j,:));
     end
     if size(timeCourse,2) > size(v,1)*32*16
         timeCourse = timeCourse(:,size(v,1)*32*16);
     end
     [S,t,f] = mtspecgramc(timeCourse',win,params);
     clear timeCourse;
     inds = f > 500;
     S = squeeze(mean(S(:,inds,:),2));
     Sf = morFilter(S',8,1250/32);
     if size(Sf,2) < size(v,1)
         Sf(:,size(Sf,2)+1:size(v,1)) = 0;
     end
%     sPlot([Sf; 10000*v(1:size(Sf,2),1).']);
%     [S,f] = mtspectrumc(S,params1);
     %size(S)
%     figure;plot(f,S);
%      S(S == 0) = nan;
%      figure;plot(f,exp(nanmean(log(squeeze(S(:,:,i))))));
%      figure;plot(f,exp(nanstd(log(squeeze(S(:,:,i))))));
end