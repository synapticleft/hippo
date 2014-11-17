
function S = getVonMisesSmoothingMatrix(pvec,mu,k)

for i=1:length(pvec)
    S(:,i) = exp(k*cos(pvec-pvec(i)+mu));
end
S = bsxfun(@rdivide,S,sum(S));