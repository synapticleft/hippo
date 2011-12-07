function odorWave(odors,data,pads,bin,data2,sniffTimes)

odors = diff(odors')';
%times = find(sniffTimes);
numOd = size(odors,1);
xDim = ceil(sqrt(numOd));
yDim = ceil(numOd/xDim);
winSize = [.5 .2];
params.Fs = 1000/bin;
params.fpass = [0 20];
params.tapers = [ceil(2/winSize(1)) winSize(1) 1];
data = [data zeros(1,pads(2))];
odors = [zeros(size(odors,1),1) odors zeros(size(odors,1),pads(2))];
figure;
for i = 1:numOd
    inds = find(odors(i,:) > 0);
    if numel(inds)
        t = [];t1 = [];
        indsNeg = find(odors(i,:) < 0);
        snippets = zeros(numel(inds),indsNeg(1)-inds(1)+range(pads)+1);
        for j = 1:numel(inds)
            %try
            snippets(j,1:(indsNeg(j)+pads(2)-inds(j)-pads(1)+1)) = data((inds(j)+pads(1)):(indsNeg(j) + pads(2)));
%            t(j) = min(times(times > inds(j)));
%            t1(j) = min(times(times > t(j)));
            %catch
            %    [i j indsNeg(j)-inds(j) size(snippets) size(data((inds(j)+pads(1)):(indsNeg(j) + pads(2))))]
            %end
        end
        if exist('data2','var')
        snippets2 = zeros(numel(inds),indsNeg(1)-inds(1)+range(pads)+1);
        for j = 1:numel(inds)
            snippets2(j,1:(indsNeg(j)+pads(2)-inds(j)-pads(1)+1)) = data2((inds(j)+pads(1)):(indsNeg(j) + pads(2)));
        end
        end
        %[S t f] = mtspecgramc(snippets',winSize,params);
        subplot(xDim,yDim,i);
        [~,inds] = sort(sum(snippets(:,-pads(1):(end-pads(2))).^2,2),'ascend');
        [~,inds] = sort(sum(max(0,snippets(:,-pads(1):(end-pads(2)))),2),'ascend');
        %[~,inds] = sort(t1-t,'ascend');
        snippets = snippets(inds,:);
        if exist('data2','var')
            snippets2 = snippets2(inds,:);
            im = zeros(size(snippets,1),size(snippets,2),3);
            snippets = max(0,snippets);
            im(:,:,1) = snippets/max(snippets(:));
            snippets2 = max(0,snippets2);
            im(:,:,2) = snippets2/max(snippets2(:));
            im(:,:,3) = im(:,:,1);
%             if i == 7 figure;
                 image(im);
%                 break
%             end
        else
        
        %figure;
            imagesc(snippets);%(t,f,squeeze(mean(S,3))');
            
        end
        allMeans(i,1:size(snippets,2)) = mean(snippets);      
        allStd(i,1:size(snippets,2)) = std(snippets);
    end
end
%sPlot(allMeans);
%sPlot(allStd);