function [A,tes,td] = allHighFreq(file,shanks,pos,all)

if all
    dec = 1;
    elecs = bsxfun(@plus,8*(ones(8,1))*(shanks-1),(1:8)');
    X = getData(file,elecs(:)',[],[],dec);
    X = filtHigh(X,1250/dec,100,8);
    %X = filtLow(X,1250,200);
    [A,tes,td] = runHighFreq(X,pos);
else
    A = [];tes = [];
    for i = shanks
        X = getData(file,(i-1)*8+(1:8));
        X = filtHigh(X,1250,100,8);
        [At,test] = runHighFreq(X,pos);
        A(:,end+1:end+size(At,2)) = At;
        tes(end+1:end+size(test,1),:,:) = test;
    end
end
h1 = figure;h2 = figure;
for i = 1:size(tes,1)
figure(h1);subplot(ceil(size(tes,1)/8),8,i);
imagesc(imfilter(squeeze(tes(i,:,:)),fspecial('gaussian',5,1),'replicate'));
axis off;
if all
    figure(h2);subplot(ceil(size(tes,1)/8),8,i);
    imagesc(reshape(sign(mean(A(:,i)))*A(:,i),[8 numel(A(:,i))/8]));axis off;
end
end
%figure;imagesc(bsxfun(@times,A,sign(skewness(A))));axis off;
if ~all
figure(h2);imagesc(bsxfun(@times,A,sign(mean(A))./max(abs(A))),[-.1 1]);axis off;
end