%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Reference
%Mike Novey and T. Adali, "ADAPTABLE NONLINEARITY FOR COMPLEX 
%      MAXIMIZATION OF NONGAUSSIANITY AND A FIXED-POINT ALGORITHM" in
% IEEE MLSP 2006.,
% Adaptable ICA algorithm based on complex generalized Gaussian distribution
% where X1 is input vector of mixtures
% and typeStr specifies MLE algorithm:
% mle_noncirc uses MLE that has noncircular model (slow)
% mle_circ assumes circular and runs mutch faster
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,W1,Z] = ACMNsym(Xin,typeStr)
inVal = 2.5;
 type = 2;
 if strcmp(typeStr,'mle_noncirc') == 1
     type = 1;
 elseif strcmp(typeStr,'mle_circ') == 1 %Much faster
     type = 2;
 end
 
 insIndex = 1;
 tol = 1e-5;
 eps = 0.001; % epsilon in G
 maxcounter = 100;

% Whitening of s:
%yyy = zeros(1,m);
% [Ex, Dx] = eig(cov(Xin'));
%   %Q = sqrt(inv(Dx)) * Ex';
% Q = mtimes(sqrt(inv(Dx)),Ex');
% 
% x = Q * Xin;

[Ex, Dx] = eig(cov(Xin'));
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
x = whiteningMatrix * Xin;
dewhiteningMatrix = Ex * sqrt (Dx);
%noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
%D = diag(flipud(1./noise_factors));
%zerophaseMatrix = E*inv (sqrt (D))*E';
rolloff_ind = 2;
noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
%D = diag(flipud(1./noise_factors));
zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';%inv (sqrt (D))*E';

  [n,m] = size(x);
 pC = (x*transpose(x))/m;

alphas1 = ones(n,1)*inVal;
%alphas1 = alphas;

% FIXED POINT ALGORITHM
  %W = zeros(n,n) + j*zeros(n,n);
%  W = complex(randn(n),randn(n));
%  Wold = zeros(n);
    W1 = complex(randn(n),randn(n));%W;
  Wold1 = zeros(n);
  k=0;
a2 = .01;
   while (norm(abs(Wold1'*W1)-eye(n),'fro')>(n*1e-4) && k < maxcounter)%&& k < 15*n)
       minAbsCos = min(abs(diag(W1' * Wold1)));
    meanAbsCos = mean(abs(diag(W1' * Wold1)));
    fprintf('Step no. %d, change in value of estimate: %.3g %.3g %.3g \n',k, norm(abs(Wold1'*W1)-eye(n),'fro'),1-minAbsCos,1-meanAbsCos);
         k = k+1;
         Wold1 = W1;
%         Wold = W;
         %%MY VERSION
         yy = W1'*x;
         localAlpha = alphas1*.5;
         absy = abs(yy).^2;
         u = absy + a2;
         u1 = bsxfun(@times,bsxfun(@power,u,localAlpha-1),localAlpha);
         %u2 = bsxfun(@times,bsxfun(@power,u,localAlpha-2),localAlpha.*(localAlpha-1));
         u2 = bsxfun(@times,u1./u,localAlpha-1);
         gRad =  x*(u1.*conj(yy)).'/size(x,2); %%double-check
         ggg = mean(u2.*absy + u1,2);
         B1 = pC*bsxfun(@times,conj(Wold1),mean(u2 .* conj(yy).^2,2).');%mean(u2.*(conj(yy).^2),2);%
         W1 = bsxfun(@times,Wold1,ggg.')-gRad + B1;%B.'*conj(Wold);
         p = alphas1;
         absy = abs(yy);
         u = absy+a2;
         up = bsxfun(@power,u,p);
         u = log(u);
         sigP = mean(bsxfun(@power,absy,p),2).^(1./p);
         gp = -(p.^(-2)).*log(p) + p.^(-2) - psi(1+1./p)./(p.^2) + ...
             mean(bsxfun(@times,up.*bsxfun(@minus,u,1./p+log(sigP)),1./((sigP.^p).*p)),2);
         ggp = 2*p.^(-3).*log(p) - 3./(p.^3) + psi(1,1+1./p).*p.^(-4)+2*psi(1+1./p).*p.^(-3) + ...
             mean(bsxfun(@times,up.*bsxfun(@plus,(u.^2 - bsxfun(@times,u,2./p)-bsxfun(@times,u,2*log(sigP))),...
             2*(p.^-2) + 2*log(sigP)./p + log(sigP).^2),1./(p.*sigP.^p)),2);
         p = p - gp./ggp;
         alphas1 = min(max(p,.2),4);
% alphas1(:) = 2.5;
%      for kk=1:n %sources
%       yy = W(:,kk)'*x; 
%       %%Fixed point
%       if type == 1  %fullAdp_mle
%           localAlpha = alphas(kk); %Because |y|^p, not (yy*)^p
%            yaug = [transpose(yy) yy']';
%           absy = abs(yy).^2;
%           u = (absy + a2);
%           u1 = localAlpha*u.^(localAlpha-1);
%           u2 = localAlpha*(localAlpha-1)*u.^(localAlpha-2);
%           gRad(:,kk) = mean(ones(n,1)*(u1.*conj(yy)).*x,2);
%           ggg = mean(u2.*absy + u1);
%           B = mean(u2.* (conj(yy).^2))*pC;
%           W(:,kk) =  Wold(:,kk)*ggg -(gRad) + (B*conj(Wold(:,kk))); 
%           [alpha ] = estimateGGDCovShapeIn(yaug, alphas(kk));
%          alphas(kk) = alpha;
%       elseif type == 2  %'fullAdp_p'
%           localAlpha = alphas(kk)*.5; %Because |y|^p, not (yy*)^p
%          absy = abs(yy).^2;
%         u = (absy + a2);
%         u1 = localAlpha*u.^(localAlpha-1);
%         u2 = localAlpha*(localAlpha-1)*u.^(localAlpha-2);
%          gRad = mean(ones(n,1)*(u1.*conj(yy)).*x,2);%/(alpha*sigP);
%          ggg = mean(u2.*absy + u1);%/(alpha*sigP);
%          B = (mean(u2.* (conj(yy).^2))*pC);%/(alpha*sigP);
%          W(:,kk) =  (Wold(:,kk)*ggg ) -(gRad) + B*conj(Wold(:,kk)); 
%          %%%Newton estimate of p
%          %yy = randn(1,1000) + j*randn(1,1000);
%          %yy = yy/std(yy);
%          alpha = alphas(kk);
%          absy = abs(yy);
%          u = (absy + a2);
%          up = u.^alpha;
%          sigP = mean(absy.^alpha)^(1/alpha); %Est nuiscance
%          p = alpha;
%          gp = -inv(p^2)*log(p) +inv(p^2) - psi(1+inv(p))/p^2 + ...
%              mean(inv(sigP^p*p)*up.*(log(u) - (inv(p) + log(sigP))*ones(1,m) ));
% %          C = mean(u.^(p-1).*log(u).*conj(yy)*inv(sigP));
%          ggp = 2*inv(p^3)*log(p) - 3*inv(p^3) + psi(1,1+inv(p))*inv(p^4)+2*psi(1+inv(p))*inv(p^3) + ...
%             mean(inv(p*sigP^p)*up.*(log(u).^2 - 2*inv(p)*log(u) -2*log(sigP)*log(u) + ...
%              (2*inv(p^2) + 2*log(sigP)*inv(p) + log(sigP)^2)*ones(1,m) ));
%          alpha = alpha - gp/ggp;
%          alpha = max(alpha,.2);
%          alpha = min(alpha,3); %2 = Gaus
%          alphas(kk) = alpha;
%       end
%   end %Loop thru sources
%     [E,D] = eig(W'*W);
%     W = W * E * inv(sqrt(D)) * E';
    [E,D] = eig(W1'*W1);
    W1 = W1 * E * inv(sqrt(D)) * E';
  end; %While
  figure;plot(alphas1);
  [mean(alphas1) median(alphas1)]
A = dewhiteningMatrix*W1;
Z = zerophaseMatrix*A;
W1 = W1' * whiteningMatrix;
%shat = W'*x;
%Ahat = inv(Q)*W;

function [alphaNew] = estimateAlphaMLEIn(yyyIn,alphaInit)
  
      yyy = abs(yyyIn);
      
      maxcounter = 10;
      tol = 1e-3;
      %alphaNew = alphaInit;
      counter = 0;
      meanY = 0;
      %%%%%%%%%%%%%%%Newton method for bivariate version
      m= size(yyy,2);
      p = alphaInit*2;
      a2 = .1;
      pNew = p;
      maxcounter = 10;
      tol = 1e-3;
      counter = 0;
      absy = (yyy).^2;
    %  while  (maxcounter - counter) >= 0;
          pNew = p;
          sigP = mean(absy.^p).^(1/p); %Est nuiscance
          %sigP = (inv(size(absy,2)-1)*sum(absy.^p)).^(1/p);
          u = (absy + a2);
          up = u.^p;

          gp = -inv(p^2)*log(p) +inv(p^2) - psi(1+inv(p))/p^2 + ...
              mean(inv(sigP^p*p)*up.*(log(u) - (inv(p) + log(sigP))*ones(1,m) ));
          %          C = mean(u.^(p-1).*log(u).*conj(yy)*inv(sigP));
          ggp = 2*inv(p^3)*log(p) - 3*inv(p^3) + psi(1,1+inv(p))*inv(p^4)+2*psi(1+inv(p))*inv(p^3) + ...
              mean(inv(p*sigP^p)*up.*(log(u).^2 - 2*inv(p)*log(u) -2*log(sigP)*log(u) + ...
              (2*inv(p^2) + 2*log(sigP)*inv(p) + log(sigP)^2)*ones(1,m) ));
          p = p - gp/ggp;
          p = max(p,.2);
          p = min(p,2);
          counter = counter+1;
          
      %end
      alphaNew = p*2-.25; %%Calibrated to p = 2 Gaussian
      alphaNew = .5*alphaNew;  %1 is Gaussian
      
function [c,R] = estimateGGDCovShapeIn(X, alphaIn)

mom=0;
converged = 1;
N = size(X,2);
R = cov(X'); %Consistant estimator but not MLE
bestC = alphaIn; %start at Gaussian
c = bestC;
Rold = zeros(2,2);
cold = 0;
delt = 10;

    xRxC = 0;
    dirXRX = 0;
    dirXRX2 = 0;
    for n = 1:N
        temp = (X(:,n)'*inv(R)* X(:,n));
        xRxC = xRxC +   real(temp^(c));
        dirXRX = dirXRX + real(log(temp)*temp^(c));
        dirXRX2 = dirXRX2 + real(log(temp)^2*temp^(c));
    end

    c2 = gamma(2/c)/(2*gamma(1/c));
    c2p = log(c2) - inv(c)*(2*psi(2/c)-psi(1/c));
    gc = N*(inv(c) - inv(c^2)*2*psi(2/c)+inv(c^2)*2*psi(1/c))-(c2^c)*(c2p*xRxC + dirXRX);

    %%Second dir
    A = N*((4*psi(2/c)/c^3) + (4*psi(1,2/c)/c^4)-(1/c^2) - (4*psi(1/c)/c^3) - (2*psi(1,1/c)/c^4));
    %Dir c2^c
    dc2C = log(c2)*(c2^c) - c*(c2^(c-1))*(c2*2*psi(2/c)/c^2 - c2*psi(1/c)/c^2);
    dc2p= -((psi(1/c) - 2*psi(2/c))/c^2) - ((psi(1,1/c) - 4*psi(1,2/c))/c^3) - ...
        ((2*psi(2/c)/c^2) - psi(1/c)/c^2);
    B = dc2C*c2p *xRxC + c2^c*(dc2p*xRxC + c2p*dirXRX);
    C = dc2C*dirXRX + c2^c*(dirXRX2);
    ggc = A - B - C; %Testted with MAthacd
    cold = c;
    cn = c-inv(ggc)*gc;
    c = min(4,max(.05,cn)); %Newton update with no negatives
    % delt = norm(Rold-R);
    delt = abs(c-cold);
    % disp([num2str(count) ' ' num2str(c)]);

  
