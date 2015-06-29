function [ stim,allIms,allPup] = eyeRegress
%find best-fit mapping of pixel eye data onto stimulus
%this is a bad idea because each movie has different coordinates!!
%but for some reason does decent!

maxLag = 30;

matFiles = {'cam2_7-9Dec16.mat','cam2_10-12Dec16.mat','cam2_13-15Dec16.mat'};
sess = [4:12];
sess1 = {4:6,7:9,10:12};
datFile = 'data_s1_corr2.mat';

load(datFile,'data','all_cam2_ons','all_cam2_offs');

sess = find(ismember([data{2:end,2}],sess));
stim = zeros(numel(sess),30);
for i = 1:numel(sess)
    inds = find(~isnan(data{sess(i)+1,40}));
    stim(i,1:min(numel(inds),size(stim,2))) = data{sess(i)+1,40}(inds(1:min(numel(inds),size(stim,2))));
end
stim = stim-stim(1,1);

allIms = zeros(numel(sess),maxLag+size(stim,2),30*40);
allPup = zeros(numel(sess),maxLag+size(stim,2));
ctr = 1;
for i = 1:numel(matFiles)
    load(matFiles{i},'pupils','ims');
    for j = sess1{i}
        for k = find(j == [data{2:end,2}])
            allIms(ctr,:,:) = ims(all_cam2_ons(k)+(0:size(allIms,2)-1),:);
            allPup(ctr,:) = data{k+1,36}(60+(0:size(allIms,2)-1));
            ctr = ctr + 1;
        end
    end
end