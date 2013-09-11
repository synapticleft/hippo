function [A,W,Z] = cfpa2(x)
%%a particular implementation of complex-ICA
%[v,whiteningMatrix,dewhiteningMatrix,zerophaseMatrix] = whiten(x);
%v = v.';
x = x.';
[N,m] = size(x);
Rxx = (x'*x)/N;
[Q,Lam] = eig(Rxx);
whiteningMatrix = Q*diag(real(diag(Lam)).^(-1/2));
dewhiteningMatrix = Q * sqrt(diag(real(diag(Lam))));
noise_factors = ones(1,size(Lam,1));
noise_factors(1:2) = .5*(1+cos(linspace(pi,0,2))); 
zerophaseMatrix = conj(Q)*sqrt(diag(fliplr(noise_factors)))*conj(Q)';
v = x*whiteningMatrix;
Phat = (transpose(v)*v)/N;
W = eye(m); 
D = W; 
%y = zeros(N,m); 
Wold = zeros(m); 
k = 0;
maxcounter = 200;

while (norm(abs(Wold'*W)-eye(m),'fro') >(m*1e-4)) && k < maxcounter %(k<15*m)
       minAbsCos = min(abs(diag(W' * Wold)));
    meanAbsCos = mean(abs(diag(W' * Wold)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g %.3g \n',k, norm(abs(Wold'*W)-eye(m),'fro'),1-minAbsCos,1-meanAbsCos);
    k = k+1;
    Wold = W;
    y = v*W;
    PhatW = Phat*W;
    for n = 1:m
        D(n,n) = transpose(W(:,n))*PhatW(:,n);
    end
    W = (v'*(y.* abs(y).^2))/N - 2*W - conj(PhatW)*D;
    [Q,Lam] = eig(W'*W);
    W = W*(Q*diag(diag(real(Lam)).^(-1/2))*Q');
    imagesc(showGrid(dewhiteningMatrix*W,[8 4]));drawnow;
end
W = (whiteningMatrix*W).';
A = pinv(W);
Z = zerophaseMatrix*A;
%A = pinv(W);
%Z = A;
