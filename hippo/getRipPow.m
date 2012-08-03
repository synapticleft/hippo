function allS = getRipPow(file,elecs,v)

d = ['/media/work/hippocampus/'];%['/media/Expansion Drive/KenjiMizuseki/'];%
file = [file '.h5'];
params.Fs = 1250;
win = [32/1250 32/1250];
params.tapers = [1/win(1) win(1) 1];params.fpass = [200 500];%[500 3000];
params1.Fs = 1250/32;params1.tapers = [15 29];
allS = zeros(numel(elecs),size(v,1));counter = 0;
info = hdf5info([d file]);
nSamples = info.GroupHierarchy.Datasets(1).Dims;
sz = nSamples(2);
figure;
for i = elecs%shanks
     X = double(h5varget([d file],'/hReal',[i-1 0],[1 sz]));
     S = log(mean((mtspecgramc(X,win,params)),2));
%      [S1,f1] = mtspectrumc(log(S),params1);
%      figure;imagesc(f,t,log(S));
%       sPlot(S',[],[],1);
%       figure;imagesc(f1,f1,log(S1));
%       figure;plot(f1,S1);
%       sPlot(S1(30:end,:)',f1(30:end),[],1);return
     Sf = morFilter(S',7.5,1250/32);
     [temp f] = mtspectrumc(S,params1);
     plot(f,log(temp));title(num2str(i));
%      hold all;
%      plot(f,log(mtspectrumc(Sf,params1)));
     drawnow;
     if numel(Sf) < size(v,1)
         Sf(numel(Sf)+1:size(v,1)) = 0;
     elseif numel(Sf) > size(v,1)
         Sf(size(v,1)+1:end) = [];
     end
     allS(counter+1,:) = Sf;
     counter = counter + 1;
end