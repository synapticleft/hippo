function [A,tes] = allHighFreq(file,shanks,pos)

A = [];tes = [];
for i = shanks
    X = getData(file,(i-1)*8+(1:8));
    X = filtHigh(X,1250,100,8);
    [At,test] = runHighFreq(X,pos);
    A(:,end+1:end+size(At,2)) = At;
    tes(end+1:end+size(test,1),:,:) = test;
end
figure;
for i = 1:size(tes,1)
subplot(size(tes,1)/8,8,i);imagesc(imfilter(squeeze(tes(i,:,:)),fspecial('gaussian',5,1),'replicate'));
axis off;
end
%figure;imagesc(bsxfun(@times,A,sign(skewness(A))));axis off;
figure;imagesc(bsxfun(@times,A,sign(mean(A))./max(abs(A))),[-.1 1]);axis off;