function shaiAnly(file,numRois)
%%INPUTS
%% file = the .lsm file (make sure you specify the directory, or set it
%%      as your 'current folder' in the matlab window).
%% numRois = the number of ROIs to select.
%% CAVEATS:
%% 1) Spine swelling causes false increase in calcium signal, therefore apparent ca 
%% of ROI is normalized by its volume (sort of like dF/F). Unfortunately this
%% introduces a slight anticorrelation between Ca and vol, which is picked up by
%% the transfer function.
%% 2) Bleaching introduces a problem because ca signal cannot predict bleach effects.
%% Current workaround is to apply a high-pass filter to volume signal to remove bleaching, but this 
%% is imperfect. Ideal would be to use a protocol that minimizes bleaching.
%% 3) The fit currently isn't cross-validated, which is prone somewhat
%% to over-fitting the data (hence the initial dip).
%% 4) ROIs arent adjusted for drift

Fs = 1.1; %sampling rate
chunkSize = 200; %number of frames to retrieve at once, if too high memory problems, if too low slow
s = tiffread302(file,11); %file that retrieves image data from .lsm
numFrames = s(1).lsm.DimensionTime;
useInds{2} = sum(s(1).data{1}) > 0;
useInds{1} = sum(s(1).data{1},2) > 0;
tempIm = zeros(sum(useInds{1}),sum(useInds{2}),3); %get rid of black regions
im = imfilter(s(1).data{2}(useInds{1},useInds{2}),fspecial('gaussian',5,.5)); %smooth
figure;imagesc((im));
axis image off;colormap gray;
for i = 1:numRois
    r{i} = roipoly;
end
fs = zeros(numRois,numFrames,2);
for i = 1:ceil(numFrames/chunkSize)
    s = getChunk(file,i,chunkSize);
    for j = 1:numel(s)
        for k = 1:2
            temp = double(s(j).data{k}(useInds{1},useInds{2}));
            for l = 1:numRois
                fs(l,(i-1)*chunkSize+j,k) = mean(temp(r{l}));
            end
           temp = imfilter(temp,fspecial('gaussian',5,1));
           tempIm(:,:,3-k) = temp/max(temp(:));
        end
        imagesc(tempIm);
        axis image off;
        drawnow;
    end
end
% for i = 1:2
%     for j =1:size(fs,1)
%         fs(j,:,i) = fs(j,:,i)/prctile(squeeze(fs(j,:,i)),20);
%     end
% end
ca = squeeze(fs(:,:,1));%./squeeze(fs(:,:,2));
vol = squeeze(fs(:,:,2));
ca = ca./vol;
subplot(211);sPlot(vol,[],0);title('Original');
vol = filtHigh(vol,Fs,.001);
subplot(212);sPlot(vol,[],0);title('Bleach-corrected (high-pass filter .001 Hz)');
%vol = bsxfun(@rdivide,vol,vol(end,:));
hist = 200;
figure;hold all;
title('transfer function');
%% xfer calculates the regression between ca and vol. 
for i = 1:size(fs,1)
    cai = zscore(ca(i,:));voli = zscore(vol(i,:));
    caT = toeplitz(zeros(hist,1),cai);
    caT = caT(:,hist+1:end);voli = voli(hist+1:end);
    xfer = (caT*caT'/size(caT,2) + eye(size(caT,1))*2)\caT*voli'/size(caT,2);
    y(i,:) = voli;
    yHat(i,:) = xfer'*caT;
    plot(xfer,'linewidth',2);hold all;
end
sPlot(complex(zscore(ca,0,2),zscore(vol,0,2)));
title('calcium (blue) and volume (red) of ROIs');
sPlot(complex(yHat,y));title('Actual (red) and predicted (blue)');

function s = getChunk(file,start,len)
s = tiffread302(file,(((start-1)*len+1):start*len)*2-1);