function X = load_datachunk(m,p,X)

rind = ceil(rand*(size(X,2) - p.imszt));
X = X(:,rind+(1:p.imszt));