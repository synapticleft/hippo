function odorClass(odors,data,bin,interval,sniffTimes,numCross,ridge,sniff,pads)
data = bsxfun(@rdivide,data,std(data,0,2));
odors = diff(odors')';
% num = 3;fasTime = 250/bin;
% f = find(sniffTimes);
% for i = 1:(numel(f)-num)
% if f(i+num)-f(i) < fasTime*(num)
% newIdx(f(i+1):f(i+num)) = 1;
% end
% end
times = find(sniffTimes);
numOd = size(odors,1);
odors = [zeros(numOd,1) odors];
allSnips = zeros(sum(odors(:)>0),size(data,1)*(range(interval)+1));
allSniffs = zeros(sum(odors(:)>0),1+range(pads));
allClasses = zeros(1,sum(odors(:)>0));
allLatent = allClasses;
counter = 0;
for i = 1:numOd
    inds = find(odors(i,:) > 0);
    if numel(inds)
        indsNeg = find(odors(i,:) < 0);
        snippets = zeros(numel(inds),size(data,1)*(range(interval)+1));
        sniffs = zeros(numel(inds),1+range(pads));
        latent = zeros(numel(inds),1);
        for j = 1:numel(inds)
            %try
            t = min(times(times > inds(j)));
            latent(j) = min(times(times > t))-t;
            %t2 = min(times(times > t1));
            t1 = t + range(interval);
            temp = zeros(size(data,1),range(interval)+1);
            temp(:,1+(0:min(range(interval),t1-t))) = data(:,t+(0:min(range(interval),t1-t)));
            te1 = temp';
            snippets(j,:) = te1(:);%min(size(snippets,2),t2-t1+1)) = data((inds(j)+pads(1)):(indsNeg(j) + pads(2)));
            sniffs(j,:) = sniff(t+(pads(1):pads(2)));
        end
        %[~,inds] = sort(sum(snippets(:,-pads(1):(end-pads(2))).^2,2),'ascend');
        %[~,inds] = sort(sum(max(0,snippets(:,-pads(1):(end-pads(2)))),2),'ascend');
        %snippets = snippets(inds,:);
    allSnips(counter+(1:numel(inds)),:) = snippets;
    allSniffs(counter+(1:numel(inds)),:) = sniffs;
    allClasses(counter+(1:numel(inds))) = i;
    allLatent(counter+(1:numel(inds))) = latent;
    counter = counter + numel(inds);
    end
end
allClasses = cumsum(min(1,diff([0 allClasses])));
figure;imagesc(allSnips);
batches = ceil(rand(size(allSnips,1),1)*3);
for i = 1:numCross
    train = batches ~= i;
    test = ~train; if numCross == 1 test = train; end
W = LDA(allSnips(train,:),allClasses(train),ridge);
L = [ones(sum(test),1) allSnips(test,:)] * W';
P = exp(L) ./ repmat(sum(exp(L),2),[1 size(L,2)]);
[~,inds] = max(P');
eff(i) = sum(inds == allClasses(test))/numel(allClasses(test));
temp = allClasses(test);
for j = 1:numel(temp)
    temp(j) = P(j,temp(j));
end
eff1(i) = mean(temp);%P(:,allClasses(~train)));
subplot(numCross,1,i);plot(inds);hold all;plot(allClasses(test));%plot(P);%
Ws(i,:,:) = W;
end
W = squeeze(mean(Ws));
L = [ones(size(allSnips,1),1) allSnips]* W';
P = exp(L) ./ repmat(sum(exp(L),2),[1 size(L,2)]);
for j = 1:numel(allClasses)
    temp(j) = P(j,allClasses(j));
end
[~,idx] = sort(temp,'descend');
figure;scatter(temp,10*allLatent);
[h x] = hist(diff(times),0:.2:50);
figure;plot(x*10,h);
figure;plot(temp);hold all;plot(temp(idx));
allSniffs = allSniffs(idx,:);
figure;plot(W(:,2:end)');
figure;imagesc(allSniffs);corr(temp',allLatent')
temp= zeros(size(allSniffs,1),size(allSniffs,2),3);
temp(:,:,1) = max(0,allSniffs);temp(:,:,3) = temp(:,:,1);
figure;image(temp/max(temp(:)));return;figure;
allClasses = allClasses(idx);
figure;
for i = 1:max(allClasses)
    subplot(3,3,i);imagesc(allSniffs(allClasses == i,:));
end
[mean(eff) mean(eff1) ]