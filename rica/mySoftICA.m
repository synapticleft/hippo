function W = mySoftICA(X)
%% Note: This sample code requires minFunc to run.
%        But you can use softICACost.m with your own optimizer of choice.
%clear all;
%addpath ~/work/RNN/minFunc/   % this should point to minFunc
                              % http://www.di.ens.fr/~mschmidt/Software/minFunc.html
                              % minFunc 2009 seems to work well
%% Load and configure a training dataset
%global params;
%params.m=20000;                 % num patches
%params.patchWidth=16;           % width of a patch
params.n=size(X,1);%size(Xparams.patchWidth^2;   % dimensionality of input

[X, mean_patch, V] = preprocess(X);
%m = sqrt(sum(patches.^2) + (1e-8));
%x = bsxfunwrap(@rdivide,patches,m);

%% Run the optimization

params.lambda = 0.05;
params.numFeatures = params.n*4;%400;
params.epsilon = 1e-3;%5

%configure minFunc
options.Method = 'lbfgs';
options.MaxFunEvals = Inf;
options.MaxIter = 300;
%options.display = 'off';
%options.outputFcn = 'showBases';

% initialize with random weights
randTheta = randn(params.numFeatures,params.n)*0.01;  % 1/sqrt(params.n);
randTheta = randTheta ./ repmat(sqrt(sum(randTheta.^2,2)), 1, size(randTheta,2)); 
randTheta = randTheta(:);

% optimize
[opttheta, cost, exitflag] = minFunc( @(theta) softICACost(theta,X, params), randTheta, options);   % Use x or xw 

% display result
W = reshape(opttheta, params.numFeatures, params.n);
display_network(W');


