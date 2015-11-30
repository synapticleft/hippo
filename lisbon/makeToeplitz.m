function X = makeToeplitz(data,lags,useZ)
if useZ
    data = zscore(data,0,2);
end
X = [];
for i = 1:size(data,1)
    X = [X; toeplitz(data(i,:),zeros(lags,1))'];
end
X = circshift(X,[0 -round(lags*.5)]);
X(end+1,:) = 1;

