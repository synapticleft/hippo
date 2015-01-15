%% various attempts at extracting features from pupil video

file = '/media/Cruzer/Roberto Data/s1_video_EMG/cam2_s1_15.avi';

% initial PCA run
m = mmread(file,2000:5:(2000+99*5));
for i = 101:-1:1
temp(i,:,:) = m.frames(i).cdata(:,:,2);
end
totFrames = abs(m.nrFramesTotal);
clear m
temp = double(temp);
mu = mean(temp(:,:));
temp(:,:) = bsxfun(@minus,temp(:,:),mu);
[u,s,v] = svd(temp(:,:),'econ');
u(end,:) = [];u(:,end) = [];
v(:,end) = [];
s(:,end) = [];s(end,:) = [];

% get compressed time series
clear temp;block = 100;
for i = 1:ceil((totFrames-1)/block)
    hiFrame = mod(min(i*block,totFrames)-1,block)+1;
    m = mmread(file,(i-1)*block+(1:hiFrame));
    for j = hiFrame:-1:1
        temp(j,:,:) = squeeze(m.frames(j).cdata(:,:,2));
    end
    temp1((i-1)*block+(1:hiFrame),:) = bsxfun(@minus,double(temp(1:hiFrame,:)),mu)*v;
    i
end

%ica
[A,W] = fastica(s*temp1,'approach','symm');

%regression
%W = pupil(pupil(:,1) ~= 0,:)'/pca100(:,pupil(:,1) ~= 0);
%W = pupil(~isnan(pupil(:,1)),:)'/pca100(:,~isnan(pupil(:,1)));

pupilOpen = pupil(:,~isnan(pupil(1,:)));
pupilOpen = bsxfun(@minus,pupilOpen,mean(pupilOpen,2));
%figure;plot(pupilOpen');
W = pupilOpen/pca100(:,~isnan(pupil(1,:)));
pupilHat = W*pca100;

%% interpolate
for i = 1:4
pupilIn(i,:) = interp1(find(~isnan(pupil(i,:))),pupil(i,~isnan(pupil(i,:))),1:size(pupil,2));
end
 pupilIn = bsxfun(@minus,pupilIn,nanmean(pupilIn,2));
WIn = pupilIn(:,~isnan(pupilIn(1,:)))/pca100(:,~isnan(pupilIn(i,:)));