function shaiAnly(file,rad)
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
chunkSize = 50; %number of frames to retrieve at once, if too high memory problems, if too low slow
s = getChunk(file,1,10);%tiffread302(file,1:10); %file that retrieves image data from .lsm
numFrames = s(1).lsm.DimensionTime;
useInds{2} = sum(s(1).data{1}) > 0;
useInds{1} = sum(s(1).data{1},2) > 0;
%tempIm(:,:,1) = imfilter(double(s(1).data{2}(useInds{1},useInds{2})),fspecial('gaussian',5,0.5)); %smooth
%tempIm(:,:,2) = imfilter(double(s(1).data{1}(useInds{1},useInds{2})),fspecial('gaussian',5,0.5)); %smooth
%tempIm(:,:,3) = 0;
%tempIm = tempIm/2^16;%max(tempIm(:));
%lastIm = squeeze(tempIm(:,:,1));mem = .5;
%lastIm = imfilter(lastIm,fspecial('gaussian',5,1));
lastIm = 0;mem = .9;
for i = 1:numel(s)
    lastIm = lastIm + double(s(i).data{2}(useInds{1},useInds{2}));
end
lastIm = lastIm/numel(s)/2^16;
fr = getFrame(size(lastIm));
rid = getRidge(lastIm,fr);
gt = graythresh(lastIm);
for i = 1:20
    rid1 = bwmorph(rid,'dilate',1);
    m = mean(lastIm(xor(rid1,rid)));
    rid = rid1;
    if m < gt*1.5
        dil = i;
        break
    end
end
%figure;imagesc(double(watershed(lastIm,1,.1)));
figure;imagesc(lastIm);%removeDen(lastIm,rid,dil));
axis image;colormap gray;
%keepGoing = 1;i = 0;
[xs ys] = ginput;
fil = fspecial('gaussian',size(lastIm),rad);
%while keepGoing %for i = 1:numRois
for i = 1:numel(xs)
%    i = i + 1;
%    temp = roipoly;
%    keepGoing = sum(temp(:));
%    if keepGoing
%        r{i} = temp;
%        temp = temp.*lastIm;
        %[~,ind] = max(temp(:));
        %ind2sub(size(temp),ind);
%        temp = regionprops(r{i},'Centroid');
        cents(i,1,:) = localMax(lastIm,fil,xs(i),ys(i));%%temp.Centroid;
end
%% junk
%cents = round(cents);
%cents(:,1,:) = [min(max(rad+1,cents(:,1,1)),size(tempIm,1)-rad) ...
%                min(max(rad+1,cents(:,1,2)),size(tempIm,2)-rad)];
%fs = zeros(numRois,numFrames,2);
% imagesc(temp);colormap gray;
% bwridgecenter(lastIm, 'bright');
% figure;
% meth = {'sobel','prewitt','roberts','log','zerocross','canny'};
% for i = 1:6
% subplot(2,3,i);
% imagesc(edge(lastIm,meth{i}));axis image;colormap gray;
%end
rid = zeros(size(lastIm));
for i = 1:ceil(numFrames/chunkSize)
    s = getChunk(file,i,chunkSize);
    for j = 1:numel(s)
        ind = (i-1)*chunkSize+j;
        if ind > 1
        for k = 1:2
            temp = double(s(j).data{k}(useInds{1},useInds{2}));
            %for l = 1:numel(r)
            %    fs(l,ind,k) = mean(temp(r{l}));
            %end
           tempIm(:,:,3-k) = temp;%imfilter(temp,fspecial('gaussian',5,1));
           %tempIm(:,:,3-k) = temp/2^16;%max(temp(:)); %%NEED TO CHANGE
        end
        tempIm(:,:,3) = 0;
        tempIm = tempIm/2^16;
        rid = getRidge(lastIm,fr);
        f = findshift(squeeze(tempIm(:,:,1)),lastIm);%'iter'
 %       tempIm = min(max(0,double(shift(tempIm,[f; 0]))),1);
        lastIm = mem*lastIm + (1-mem)*squeeze(tempIm(:,:,1));
        lastImDen = removeDen(lastIm,rid,dil+1);
        subplot(121);imagesc(lastImDen);hold all;axis image;colormap gray;
        subplot(122);imagesc(rid);axis image;colormap gray;
        %subplot(133);imagesc(tempIm);
        subplot(121);
         for k = 1:size(cents,1)
             cents(k,ind,:) = localMax(lastImDen,fil,cents(k,ind-1,1),cents(k,ind-1,2));
%             sh = round(squeeze(cents(k,ind-1,:))'-fliplr((size(r{1})+1)/2));
%             temp = lastIm.*circshift(fil,[sh(2) sh(1)]);%squeeze(tempIm(:,:,1))
%             [~,maxInd] = max(temp(:));
%             [cents(k,ind,2) cents(k,ind,1)] = ind2sub(size(temp),maxInd);
% %            temp = cents(k,ind-1,:);
% %            bounds = [temp(1)-rad temp(1)+rad temp(2)-rad temp(2)+rad];%[max(1,temp(1)-rad) min(size(tempIm,1),temp(1)+rad) ...
% %                %max(1,temp(2)-rad) min(size(tempIm,2),temp(2)+rad)];
% %            f = findshift(squeeze(tempIm(bounds(1):bounds(2),bounds(3):bounds(4),1)),...
% %                squeeze(lastIm(bounds(1):bounds(2),bounds(3):bounds(4),1)),'integer');
% %            cents(k,ind,:) = [min(max(rad+1,cents(k,1,1)+f(2)),size(tempIm,1)-rad) ...ind-1
% %                min(max(rad+1,cents(k,1,2)+f(1)),size(tempIm,2)-rad)];%ind-1
             scatter(cents(k,ind,1),cents(k,ind,2),'filled');
        end
        hold off;
        drawnow;
        end
    end
end
return
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
his = 200;
figure;hold all;
title('transfer function');
%% xfer calculates the regression between ca and vol. 
for i = 1:size(fs,1)
    cai = zscore(ca(i,:));voli = zscore(vol(i,:));
    caT = toeplitz([cai(1); zeros(his-1,1)],cai);
    caT = caT(:,his+1:end);voli = voli(his+1:end);
    xfer = (caT*caT'/size(caT,2) + eye(size(caT,1))*2)\caT*voli'/size(caT,2);
    y(i,:) = voli;
    yHat(i,:) = xfer'*caT;
    plot(xfer,'linewidth',2);hold all;
end
sPlot(complex(zscore(ca,0,2),zscore(vol,0,2)));
title('calcium (blue) and volume (red) of ROIs');
sPlot(complex(yHat,y));title('Actual (red) and predicted (blue)');

function f = getFrame(sz)
f = zeros(sz);
f([1 sz(1)],:) = 1;
f(:,[1 sz(2)]) = 1;

function xy = localMax(im,fil,x,y)
    sh = round([x y] - (size(im')+1)/2);
    %        sh = round(squeeze(cents(k,ind-1,:))'-fliplr((size(r{1})+1)/2));
    %subplot(211);imagesc(im);colormap('gray');
    im = im.*circshift(fil,[sh(2) sh(1)]);
    %subplot(212);
    %imagesc(im);pause(.5);colormap('gray');
    [~,maxInd] = max(im(:));
    [xy(2) xy(1)] = ind2sub(size(im),maxInd);

function lastIm = removeDen(lastIm,rid,dil)
rid = imfilter(double(bwmorph(rid,'dilate',dil)),fspecial('gaussian',5,2),'replicate');
lastIm = lastIm.*(1-rid);

function im = getRidge(im,fr)
im = xor(fr,double(watershed(imfilter(im,fspecial('gaussian',5,1)),1,.05)));%
% temp = bwmorph(im2bw(im,graythresh(im)),'close','inf');
% temp = bwmorph(temp,'thin','inf');
% while 1%s > 2%i = 1:10
%     e = bwmorph(temp,'endpoints');
%     if sum(e(:)) > 2
%         %s = sort(im(e),'descend');
%         %e(e) = im(e) < s(2);
%         temp = temp - e;
%     else
%         break;
%     end
% end
% e = bwmorph(temp,'branchpoints');
%if sum(e(:))
%    temp = rid;
%end

function c = getCentroids(s)
for k = 1:numel(s)
    idx = s(k).PixelIdxList;
    pixel_values = double(I(idx));
    sum_pixel_values = sum(pixel_values);
    x = s(k).PixelList(:, 1);
    y = s(k).PixelList(:, 2);
    xbar = sum(x .* pixel_values) / sum_pixel_values;
    ybar = sum(y .* pixel_values) / sum_pixel_values;
end

function s = getChunk(file,start,len)
s = tiffread302(file,(((start-1)*len+1):start*len)*2-1);