%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference
% Mike Novey and T. Adali, "On Extending the complex FastICA algorithm
% to noncircular sources" IEEE Trans. Signal Processing, vol. 56, no. 5, pp. 2148-2154, May 2008.
%
% Performs symmetric orthogonalization of the
% non-circular complex FastICA algorithm where
% xold is the mixtures and typeStr = 'log', 'kurt', or 'sqrt' is
% the nonlinearity.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [A, W1] = nonCircComplexFastICAsym(xold,typeStr)
type = 0;
if strcmp(typeStr,'log') == 1
    type = 1;
elseif strcmp(typeStr,'kurt') == 1
    type = 2;
elseif strcmp(typeStr,'sqrt') == 1
    type = 3;
end

tol = 1e-3;
a2= 1;
maxcounter = 200;
% Whitening of s:
% yyy = zeros(1,m);
% [Ex, Dx] = eig(cov(xold'));
% Q = mtimes(sqrt(inv(Dx)),Ex');
% x = Q * xold;

[Ex, Dx] = eig(cov(xold'));
d = flipud(diag(Dx));
cumVar = sum(d);
maxLastEig = sum(cumsum(d)/cumVar < .9999999)
Dx = Dx(end-maxLastEig+1:end,end-maxLastEig+1:end);
Ex = Ex(:,end-maxLastEig+1:end);
factors = diag(Dx);
noise_factors = ones(size(Dx,1),1);
rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .999999)
noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
Dx = diag(factors./noise_factors);
whiteningMatrix = sqrt(inv(Dx)) * Ex';
x = whiteningMatrix * xold;
dewhiteningMatrix = Ex * sqrt (Dx);

[n,m] = size(x);
%Pseudo-covariance
pC = (x*transpose(x))/m;
% FIXED POINT ALGORITHM
W1 = complex(randn(n),randn(n));%eye(n);%
Wold1 = zeros(size(W1));%zeros(n);
k=0;

while (norm(abs(Wold1'*W1)-eye(n),'fro')>(n*tol) && k < maxcounter)%(min(abs(diag(W1' * Wold1))) > (n*tol)
    minAbsCos = min(abs(diag(W1' * Wold1)));
    meanAbsCos = mean(abs(diag(W1' * Wold1)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g %.3g \n',k, norm(abs(Wold1'*W1)-eye(n),'fro'),1-minAbsCos,1-meanAbsCos);
    k = k+1;
    Wold1 = W1;
    %%MY VERSION
%    tic;
    yy = Wold1'*x;
    absy =abs(yy).^2;
    %%Fixed point
    if type == 1 %%log
        g = 1./(a2 + absy);
        gp =  -1./(a2 + absy).^2;
    elseif type == 2  %Kurt
        g = absy;
        gp =  ones(size(absy));
    elseif type == 3  %sqrt
        g = 1./(2*sqrt(a2 + absy));
        gp =  -1./(4*(a2 + absy).^(3/2));
    end
    gRad =  x*(g.*conj(yy)).'/size(x,2);%mean(ones(n,1)*(g.*conj(yy)).*x,2);
    ggg = mean(gp.*absy + g,2);
    B = pC*bsxfun(@times,conj(Wold1),mean(gp .* conj(yy).^2,2).');%*pC;
    W1 = bsxfun(@times,Wold1,ggg.')-gRad + B;
    [E,D] = eig(W1'*W1);
    W1 = W1 * E / (sqrt(D)) * E';
%    toc
%     %%ORIGINAL VERSION
%     tic;
%     for kk=1:n %Loop thru sources
%         yy = W(:,kk)'*x;
%         absy =abs(yy).^2;
%         %%Fixed point
%         if type == 1 %%log
%             g = 1./(a2 + absy);
%             gp =  -1./(a2 + absy).^2;
%         elseif type == 2  %Kurt
%             g = absy;
%             gp =  ones(size(absy));
%         elseif type == 3  %sqrt
%             g = 1./(2*sqrt(a2 + absy));
%             gp =  -1./(4*(a2 + absy).^(3/2));
%         end
%         gRad =  mean(ones(n,1)*(g.*conj(yy)).*x,2);
%         ggg = mean(gp.*absy + g);
%         B = mean(gp .* conj(yy).^2)*pC;
%         W(:,kk) =  Wold(:,kk)*ggg -(gRad) + (B*conj(Wold(:,kk)));
%     end
%     %Orthonormalization
%     [E,D] = eig(W'*W);
%     W = W * E * inv(sqrt(D)) * E';
%     toc
%     1-corr(W(:),W1(:))
%     scatter(real(W(:)),real(W1(:)));hold all;scatter(imag(W(:)),imag(W1(:)),'r');drawnow;hold off;
end; %Loop thru sources

%shat = W1'*x; %Estimated sources
A = dewhiteningMatrix*W1;
W1 = W1' * whiteningMatrix;

%Ahat = inv(Q)*W1;

