function [A,W,Z] = cfpa3(x)
%a particular implementation of complex ICA -- not sure if it works right?
[m,N] = size(x);
Rxx = (x*x')/N;
[Ex,Dx] = eig(Rxx);
whiteningMatrix = sqrt(inv(Dx)) * Ex';
v = whiteningMatrix * x;
Phat = v*transpose(v)/N;

%dewhiteningMatrix = Ex * sqrt (Dx);
rolloff_ind = 2;
noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';%inv (sqrt (D))*E';

W = eye(m); 
D = W; 
Wold = zeros(m); 
k = 0;
maxcounter = 20;

while (norm(abs(Wold'*W)-eye(m),'fro') >(m*5e-4)) && k < maxcounter %(k<15*m)
       minAbsCos = min(abs(diag(W' * Wold)));
    meanAbsCos = mean(abs(diag(W' * Wold)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g %.3g \n',k, norm(abs(Wold'*W)-eye(m),'fro'),1-minAbsCos,1-meanAbsCos);
    k = k+1;
    Wold = W;
    y = W*v;
    PhatW = W*Phat;
    for n = 1:m
        D(n,n) = transpose(W(:,n))*PhatW(:,n);
    end
    W = (v'*(y.* abs(y).^2))/N - 2*W - conj(PhatW)*D;
    [Q,Lam] = eig(W'*W);
    W = W*(Q*diag(diag(real(Lam)).^(-1/2))*Q');
end
W = (whiteningMatrix*W).';
A = pinv(W);
Z = zerophaseMatrix*A;
%A = pinv(W);
%Z = A;
