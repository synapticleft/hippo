function GLMfits(file)
%% fit all cells in a run using multi-resolution wavelet, multi-electrode LFP data
%ec013 - [1560096 15 64]
%ec014 - [ 15 128]
d = '/media/work/hippocampus/ec013.670/';
dec = 32;
bin = 1/1.25;
dims = [30 1560096 15];
someShanks = 5:8;

levels = 1:.5:8;
pos = [d file '.whl'];
[~,pos,isFast,~] = fixPos(pos);
%runNum = round(interp(runNum,dec));
isFast = logical(round(interp([double(isFast); 0],dec)));
[spikeTime,spikeID,~,info] = LoadCluRes([d file],someShanks);
spikeTime = ceil(spikeTime/info.SampleRate*1000/bin);
[~, cellInfo] = hipSpikes(file,bin,d,1);
spikeTime(ismember(spikeID,find(cellInfo.type < 1))) = [];
spikeID(ismember(spikeID,find(cellInfo.type < 1))) = [];
%allIDs = unique(spikeID);
h = hist(spikeID,1:max(spikeID));
allIDs = find(h > 500);
allIDs = [3 5 28];
m = memmapfile([[d file] 'CwtSpk.dat'],'Format',{'single' dims 'X'});

y = zeros(dims(2),1);

%trainSamples = false(size(y));
%testSamples = trainSamples;
splitData = rand(1,sum(isFast)) > .5;
trainSamples = splitData;%(isFast)
testSamples= ~splitData;
%binned = zeros(numel(levels),numel(allIDs),
levelsMine = [1 2 7 13 14];%[1 1.5 4 7 7.5];
for i = 1:numel(levelsMine)
    tic;
    X = permute(m.Data(1).X([levelsMine(i) levelsMine(i)+numel(levels)],isFast,:),[2 1 3]);
    toc
    for j = 1:numel(allIDs)
        y(:) = 0;
        y(spikeTime(spikeID == allIDs(j))) = 1;
        y = y(isFast);
        tic;
        X = zscore(X(:,:));
        %gg(i,j,:) = glmfit(X(trainSamples,:),y(trainSamples),'binomial');%'normal');%
        %sim = glmval(squeeze(gg(i,j,:)),X(testSamples,:),'logit');%'identity');%%
        ggL(i,j,:) = y(trainSamples)'/X(trainSamples,:)';
        sim = X(testSamples,:)*squeeze(ggL(i,j,:));
        cc(i,j) = corr(y(testSamples),sim);
        %sim = X*squeeze(ggL(i,j,:));
        %binned(i,j,:,:) = accumarray([
        toc
        save([d 'allFitsLin.mat'],'cc','ggL','allIDs');
        load([d 'allFits.mat'],'gg');
        sims(i,j,:) = glmval(squeeze(gg(i,j,:)),X,'logit');
    end
end
save([d 'allFits.mat'],'sims','-append');


%%%%%%%%%%%%%%%%%%%%%%
% global RefreshRate;  % Stimulus refresh rate (Stim frames per second)
% RefreshRate = 1250; 
% DTsim = 1;
% gg0 = gmakeFittingStruct_GLM(zeros(taps,size(Stim,2)),DTsim);  % projects sta into basis for fitting k
% gg0.nlfun = @logexp1;
% gg0.tsp = tsp;  % Insert spikes into fitting struct
% gg0.tspi = 1;   % First spike to use (you can ask it to ignore the first "n" spikes)
% %[logli0,rr0,tt] = neglogli_GLM(gg0,Stim); % Compute logli of initial params (if desired)
% 
% % Do ML estimation of model params
% opts = {'display', 'iter', 'maxiter', 100};
% tic;
% [gg, negloglival] = MLfit_GLM(gg0,Stim,opts); % do ML (requires optimization toolbox)
% toc

% Stim1 = zeros(numel(levels)*size(Stim,1),size(Stim,2));
% tic;
% for i = 1:size(Stim,1)
%   Stim1((i-1)*numel(levels)+(1:numel(levels)),:) = cwt(Stim(i,:),2.^(levels),'cmor1-1');
% end
% if ~isreal(Stim1)
%     Stim1 = [real(Stim1);imag(Stim1)];
% end
% Stim1 = zscore(Stim1,0,2);
% toc
% for i = 1:numel(levels)%size(Stim,1)
%     Stim1((i-1)*size(Stim,1)+(1:size(Stim,1)),:) = circshift(Stim,[0 round(i-numel(levels)/2)]);%cwt(Stim(i,:),2.^(levels),'morl');
% end