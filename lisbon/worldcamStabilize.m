function allData1 = worldcamStabilize(inFile,ff)


vidReader = VideoReader(inFile);
if ~exist('ff','var')
    vidReader.currentTime = vidReader.Duration/2;
    frameRGB = double(readFrame(vidReader));
    imagesc(frameRGB/256);
    r = getrect
    ff = frameRGB(round(r(2)+(0:r(4))),round(r(1)+(0:r(3))),:);
    vidReader.currentTime = 0;
end

allData = [];
while hasFrame(vidReader)
    frameRGB = double(readFrame(vidReader));
    fff = imfilter(frameRGB(:,:,1),ff(:,:,1))-imfilter(frameRGB(:,:,2),ff(:,:,2))*3;
    lm = vision.LocalMaximaFinder('MaximumNumLocalMaxima',4,'IndexDataType','double','Threshold',500000,'NeighborhoodSize',round([size(ff,2) size(ff,1)]/2)*2+1); 
    s = step(lm,fff);
    s(size(s,1)+1:4,1:2) = NaN;
    allData = [allData; s(:)'];% step(lm,fff);
    %imagesc(frameRGB/256);hold on;
    %scatter(s(:,1),s(:,2),'k','filled');drawnow;hold off;
end
allData1 = nan*ones(size(allData));
for i = 1:size(allData,1)
    if sum(isnan(allData(i,:))) == 0
        [~,ind] = sort(allData(i,1:4)+allData(i,5:8));
        allData1(i,[1 4 5 8]) = allData(i,[ind(1) ind(4) ind(1)+4 ind(4)+4]);
        [~,ind] = sort(allData(i,1:4)-allData(i,5:8));
        allData1(i,[2 3 6 7]) = allData(i,[ind(1) ind(4) ind(1)+4 ind(4)+4]);
    end
end

csvwrite([inFile(1:end-4) 'Points.csv'],allData1);