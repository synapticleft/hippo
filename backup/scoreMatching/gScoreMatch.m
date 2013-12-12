function [A,W,Z] = gScoreMatch(Xin,nHidden)

[Ex, Dx] = eig(cov(Xin'));
d = flipud(diag(Dx));
cumVar = sum(d);
maxLastEig = numel(d);%sum(cumsum(d)/cumVar < .9999999)
Dx = Dx(end-maxLastEig+1:end,end-maxLastEig+1:end);
Ex = Ex(:,end-maxLastEig+1:end);
factors = diag(Dx);
noise_factors = ones(size(Dx,1),1);
%rolloff_ind = sum(cumsum(flipud(factors))/cumVar > .999999)
%noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi-.01,0,rolloff_ind))); 
Dx = diag(factors./noise_factors);
whiteningMatrix = sqrt(inv(Dx)) * Ex';
whiteData = whiteningMatrix * Xin;
dewhiteningMatrix = Ex * sqrt (Dx);
rolloff_ind = 2;
%noise_factors(:) = 1;
%noise_factors(1:rolloff_ind) = .5*(1+cos(linspace(pi,0,rolloff_ind)));
zerophaseMatrix = Ex*sqrt(diag(flipud(noise_factors)))*Ex';

%% Run the optimization with minFunc (ICA)
fprintf('\nTraining ICA (w/ Score Matching)\n\n');
nInput = size(whiteData, 1);
W = randn(nHidden, nInput); 
options.Method  = 'lbfgs';
options.maxIter = 100;	    % Maximum number of iterations of L-BFGS to run 
options.display = 'on';

tic
[W, cost] = minFunc( @icaScoreMatching, ...
                        W(:), options, whiteData, ...
                        nHidden, nInput);
toc

%% Display Results
W = reshape(W, nHidden, nInput);
A = dewhiteningMatrix*W';
Z = zerophaseMatrix*A;
W = W * whiteningMatrix;
