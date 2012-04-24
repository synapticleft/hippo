function allData = golgi
%normalizing by max is problematic because gfp max differs for Cu and BCS
%(check out norms(red)/norms(green).
%figure;
d = dir;
isub = [d(:).isdir]; %# returns logical vector
nameBatches = {d(isub).name}';
nameBatches(ismember(nameBatches,{'.','..'})) = [];
allData = [];
for batch = 1:numel(nameBatches)
    cd(nameBatches{batch});
    gd = dir;
    isub = [gd(:).isdir]; %# returns logical vector
    gd = {gd(isub).name}';
    gd(ismember(gd,{'.','..'})) = [];
    for gene = 1:numel(gd)
        allData(end+1).batchName = nameBatches{batch};
        allData(end).geneName = gd{gene};
        cd(gd{gene});
        d = dir('*1.TIF');
        norms = zeros(numel(d),4);
        for i = 1:numel(d)
            for j = 1:2
                x(:,:,j) = double(imread([d(i).name(1:end-5) num2str(j+1) '.TIF']));
                temp = x(:,:,j);
                h = filtfilt(gausswin(5),1,hist(temp(:),1:1001));
                h = h(1:end-1);
                [~,p] = findpeaks(h/max(h(:)),'MINPEAKHEIGHT',.5);
                try
                    norms(i,j) = p(1);
                catch
                    plot(h/sum(h(:)));
                end
                temp = temp - p(1);
                norms(i,j+2) = prctile(temp(:),98);
                temp = temp/norms(i,j+2);
                x(:,:,j) = temp;
            end
            x = x(:,:,[2 1]);
            f = findshift(squeeze(x(:,:,1)),squeeze(x(:,:,2)),'grs');
            x(:,:,2) = shift(squeeze(x(:,:,2)),f);
            nls(i,:) = nlImage(x);
            %[b(i,:) cc(i)] = glmImage(x,rad);
        end
        hold off;
        allData(end).norms = norms;
        allData(end).nls = nls;
        cd ..;
    end
    cd ..;
end

function y = nlImage(im)
im = imresize(im,1/3);
x = squeeze(im(:,:,1));x = x(:);
y = squeeze(im(:,:,2));y = y(:);
%subset = rand(numel(y),1) < .1;
%f = fit(y,x,'rat11');
%plot(f,y(subset),x(subset));pause(1);
xBin = round((x+.1)*10);xBin = max(0,min(12,xBin))+1;
x = accumarray(xBin,x,[13 1],@mean);
y = accumarray(xBin,y,[13 1],@mean);
plot(y);hold all;drawnow;

function [b cc] = glmImage(im,rad)
im = imresize(im,1/3);thresh = .15;
%try 'identity', adding blue channel, diff radii
numCross = 3;link = 'identity';%'logit';%
X = zeros((2*rad+1)^2,size(im,1)*size(im,2));
y = squeeze(im(:,:,2));y = y(:);
counter = 1;
for i = -rad:rad
    for j = -rad:rad
        temp = circshift(squeeze(im(:,:,1)),[i j]);
        X(counter,:) = temp(:);
        counter = counter + 1;
    end
end
y = y(X((size(X,1)-1)/2,:) > thresh);
X = X(:,X((size(X,1)-1)/2,:) > thresh);
rInd = ceil(rand(numel(y),1)*numCross);
for i = 1:numCross
    b(i,:) = glmfit(X(:,rInd ~= i)',y(rInd ~= i),'normal','link',link);
    temp = glmval(b(i,:)',X(:,rInd == i)',link);
    cc(i) = corr(temp,y(rInd == i));
end
b = mean(b);
subplot(221);imagesc(clipIm(im));axis image off;
%im(:,:,1) = reshape(glmval(b',X',link),[size(im,1) size(im,2)]);
%subplot(222);imagesc(clipIm(im));axis image off;
subplot(223);imagesc(real(reshape(b(2:end),(2*rad+1)*[1 1])));colorbar;axis image off;
subplot(224);imagesc(squeeze(im(:,:,2)) > thresh);axis image off;drawnow;
cc = mean(cc);

function im = clipIm(im)
im = min(1,max(0,real(im)));


%     subplot(221);imagesc(power(max(0,min(1,x)),1/2));axis image off;
%     x1 = permute(x,[3 1 2]);
%     subplot(223);imagesc(bins{1},bins{2},sqrt(hist3(x1(1:2,:)',bins)));
%     x = idealfilterG(x,1,1);%imfilter(x,fspecial('gaussian',10,4));
%     x1 = permute(x,[3 1 2]);
%     temp = squeeze((x1(3,:,:))) > .2;
%     %[u s v] = svd(bsxfun(@minus,x1(1:2,temp),mean(x1(1:2,temp),2)),'econ');
%     %binsa{1} = linspace(-3*std(v(:,1)),3*std(v(:,1)),100);binsa{2} = binsa{1};
%     %figure;imagesc(temp);return
%     %figure;imagesc(binsa{1},binsa{2},sqrt(hist3(v,binsa)));return
%     subplot(222);imagesc(power(max(0,min(1,x)),1/2));axis image off;
%     subplot(224);imagesc(bins{1},bins{2},sqrt(hist3(x1(1:2,:)',bins)));pause(.1);return