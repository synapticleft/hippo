function [im frameCol] = superImp(x,frames,rad,maxVal,numCol,bb)
%%SUPERIMP combines multiple components into 1 image, assigning each
%%component a different color, for each pixel, choosing the component with
%%the largest magnitude at that location. All components are normalized.
%%INPUTS:   c = all 2-d image components
%%          frames = which components to combine in image (choose ones with well-defined features)
%%          rad = width of gaussian smoothing kernel
%%          maxVal = if you want to normalize all components by a fixed value (default: normalize maximum of each component to 1)
%%          bb = in case you want a separate background image to assign the intensity of the pixel, rather than the componnent
% if exist('frames','var') && ~isempty(frames)
%     x = x(frames(randperm(numel(frames))),:,:);
% else
%     x = x(randperm(size(x,1)),:,:);
% end
if exist('frames','var')
    x = x(frames,:,:);
end
allReg = zeros(size(x));
for i = 1:size(x,1)
    if exist('rad','var') && rad
        x(i,:,:) = imfilter(squeeze(x(i,:,:)),fspecial('gaussian',5,rad));
    end
    allReg(i,:,:) = x(i,:,:) > max(max(squeeze(x(i,:,:))))/2;%maxVal/4;%
    if ~exist('maxVal','var')
        x(i,:,:) = x(i,:,:)/max(max(max(x(i,:,:))));
    else
        x(i,:,:) =  x(i,:,:)/maxVal;
    end
end
if exist('numCol','var') && ~isempty(numCol)
    cols = sepCols(allReg,numCol);
    cols = cols/numCol;
else
    cols = (1:size(allReg,1))/size(allReg,1);
end
x = abs(x);
x = min(1,max(0,x));
[a b]= max(x);
a = squeeze(a); b = squeeze(b);
im = zeros(3,size(a,1),size(a,2));
for i = 1:size(x,1);
    im(1,b == i) = cols(i);%i/size(x,1);
    frameCol(1,i,1) = cols(i);
    im(2,b == i) = 1;
    im(3,b == i) = a(b == i);
    if exist('bb','var')
        bb = squeeze(mean(bb));
        bb = bb-150;
        bb = bb./1500;bb = max(0,min(bb,1));%max(max(bb));
        im(3,b == i) = bb(b == i);
    end
end
frameCol(1,:,2:3) = 1;
frameCol = hsv2rgb(frameCol);
im = max(0,im);
im = permute(im,[2 3 1]);
im = hsv2rgb(im);
figure;image(im);

function s = sepCols(reg,numCol)
for i = 1:size(reg,1)
    imagesc(squeeze(reg(i,:,:)));pause(.03);
    allValid = 1:numCol;
    for j = 1:i-1
        if sum(sum(sum(reg(i,:,:).*reg(j,:,:))))
            allValid(allValid == s(j)) = [];
        end
    end
    temp = randperm(numel(allValid));
    s(i) = allValid(temp(1));
end