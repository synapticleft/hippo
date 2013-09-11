function day = allLinear(ratID)
%%This was used to pull all spiking and LFP data from Kenji's linear files,
%%and write them to respective matlab files.
writeDir = '/media/Seagate Expansion Drive/';
infoDir = '/media/Expansion Drive/redwood/';

shankInfo = importdata([infoDir 'ElePosition.txt']);
load([infoDir 'KenjiData.mat'],'Beh');
counter = 1;
for i = 1:numel(shankInfo)
    if ~isempty(strfind(shankInfo{i},ratID)) && ~isempty(strfind(shankInfo{i},'CA1'))
        validShank{counter} = shankInfo{i}(1:find(isspace(shankInfo{i}),1)-1);
        counter = counter + 1;
    end
end
day = [];
for i = 1:size(Beh,1)
if ~isempty(strfind(Beh{i,5},'linear')) && strcmp(Beh{i,8},'')
    [tf,loc] = ismember(Beh{i,1},validShank);
%     if tf
%     load([writeDir 'gaMatlab/' Beh{i,4} '.mat'],'cellInfo');
%     tf = tf & min(cellInfo.ID) < 0;
%     end
    if tf
%         fileBase = [writeDir 'kenji_data/' Beh{i,1} '/' Beh{i,4} '/' Beh{i,4}];
%         X = getData(fileBase,32);
%         Xf = morFilter(X,8,1250/32);
%         [u,~] = eig(Xf*Xf');
%         v = (u(:,end)\Xf)';
%         pos = importdata([fileBase '.whl']);
%         day = [day; [i loc]];
        [sp, cellInfo] = hipSpikes(Beh{i,4},32/1.25,[writeDir 'kenji_data/' Beh{i,1} '/' Beh{i,4} '/']);%,2);
%        spf = morFilter(sp,8,1250/32);
%        save([writeDir 'gaMatlab/' Beh{i,4} '.mat'],'X','Xf','v','pos','sp','spf','cellInfo','shankInfo');
        save([writeDir 'gaMatlab/' Beh{i,4} '.mat'],'cellInfo','-append');
        printf([Beh{i,4} '.mat\n']);
    end
end
end
% save([writeDir 'gaMatlab/linearDays.mat'], 'day'); 

%% get .whl, downsampled .eeg, spikes, filter each, 1st PC, demodulate, decoding matrix

