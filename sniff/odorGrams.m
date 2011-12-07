function odorGrams(odors,data,pads,bin)

odors = diff(odors')';

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
        indsNeg = find(odors(i,:) < 0);
        snippets = zeros(numel(inds),indsNeg(1)-inds(1)+range(pads)+1);
        for j = 1:numel(inds)
            try
            snippets(j,1:(indsNeg(j)+pads(2)-inds(j)-pads(1)+1)) = data((inds(j)+pads(1)):(indsNeg(j) + pads(2)));
            catch
                [i j indsNeg(j)-inds(j) size(snippets) size(data((inds(j)+pads(1)):(indsNeg(j) + pads(2))))]
            end
        end
        [S t f] = mtspecgramc(snippets',winSize,params);
        subplot(xDim,yDim,i);
        imagesc(t,f,squeeze(mean(S,3))');
        allMeans(i,1:size(snippets,2)) = mean(snippets);        
    end
end
sPlot(allMeans);