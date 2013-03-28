function [h nFields momInertia momInertia1] = getPlaceStats(filts,acts,fullResp,thresh)% momInertia cents

momInertia = [];
%angChange = [];
nFields = [];
cents = [];
%counter = 0;
%accumResp = zeros(size(fullResp,2),size(fullResp,3));
%for i = 1:numel(filts)
%    for j = 1:size(fullResp,2)
%        if numel(acts{i})
%            accumResp(j,squeeze(acts{i}(~isnan(acts{i}(:,j,1)),j,1))) = 1;
%        end
%    end
%end
%accumResp = mean(accumResp);
%params.tapers = [3 5];params.Fs = (size(fullResp,3)/2)/250;
%accumResp = reshape(accumResp,[numel(accumResp)/2 2]);
%accumResp = bsxfun(@minus, accumResp,mean(accumResp));
%[S,f] = mtspectrumc(accumResp,params);

for i = 1:numel(filts)
    temp = filts{i};
    nFields(i) = size(temp,2);
    for j = 1:size(temp,2)
        grid = (1:size(temp,1))';
        cm = sum(grid.*abs(temp(:,j)))/sum(abs(temp(:,j)));
        momInertia = [momInertia sum(abs((grid-cm)).*abs(temp(:,j)))/sum(abs(temp(:,j)))];
%        angChange = [angChange angle(temp(2:end,j).'/temp(1:end-1,j).')];
    end
end

range = 20;
xs = meshgrid(1:size(fullResp,3),1:size(fullResp,2));
bins{1} = -range/2:range/2;
angBins = 20;
bins{2} = linspace(-1,1,angBins)*(angBins-1)/angBins*pi;
h = 0;
allDat = [];
fullResp1 = fullResp;fullResp1(:) = zscore(fullResp(:));
xdim = 6;ydim = 7;counter = 1;
for i = 1:numel(filts)
    if numel(acts{i})% && size(filts{i},2) <= 5
        meanLocs = round(nanmean(acts{i}(:,:,1),2));
        for j= 1:numel(meanLocs)
            inds= meanLocs(j)+(-range/2:range/2);
            inds(inds <= 0) = [];inds(inds > size(fullResp,3)) = [];
            tempIm = squeeze(fullResp(i,:,inds));
            absResp = abs(squeeze(fullResp1(i,:,inds)));
            tempXs = xs(:,inds)-meanLocs(j);
            allDat = [allDat ;[tempXs(absResp > thresh) tempIm(absResp > thresh)]];
            %subplot(xdim,ydim,counter);imagesc(hist3([tempXs(absResp > thresh) angle(tempIm(absResp > thresh))],bins)');axis off;drawnow;
            counter = counter + 1;
        end
        cents = [cents; meanLocs];
    end
end
%Fs = (size(fullResp,3)/2)/250;
ac = accumarray(allDat(:,1)-min(allDat(:,1))+1,allDat(:,2),[],@mean);
%sPlot(ac.')
%allDat = bsxfun(@minus,allDat,mean(allDat));
allDat(:,2) = angle(allDat(:,2));%*exp(1i*mean(allDat(:,2))));
h = hist3(allDat,bins);
%[beta betaInt] = regress(angle(ac),bins{1}.')
%figure;imagesc(bins{1},bins{2},h');drawnow;hold all;plot(bins{1},angle(ac),'w','linewidth',2);%plot(bins{1},bins{1}*beta,'w','linewidth',2);%
%figure;plot(f,S);
%figure;hist(angChange*Fs/(2*pi));
%figure;hist(nFields,0:10);

fullResp = permute((fullResp),[1 3 2]);
fullResp = fullResp(nFields>0,:,:);
temp = circshift(fullResp,[0 0 1]);
buff = 12;%range*2;
selfCorr = zeros(size(fullResp,1),2*buff+1);
allCorr = zeros(size(fullResp,1),size(fullResp,1),2*buff+1);
momInertia1 = [];
for i = 1:size(fullResp,1)
    selfCorr(i,:) = xcov(fullResp(i,:),temp(i,:),buff,'coeff');
    %for j =1:size(fullResp,1)%i-1
    %    allCorr(i,j,:) = xcov(fullResp(i,:),fullResp(j,:),buff,'coeff');
    %end
    %allCorr(i,i,:) = 0;%selfCorr(i,:);
    %[~,m] = max(max(abs(squeeze(allCorr(i,:,buff+1:end))),[],2));
    %allCorr(i,~ismember(1:size(fullResp,1),m),:) = 0;
    momInertia1 = [momInertia1 sum(abs((grid-(size(selfCorr,2)-1)/2)).*abs(selfCorr(i,:))')/sum(abs(selfCorr(i,:)'))];
end