load('dat1_b0500')
%% Example of 1D decoding with various methods (Template Matching, OLE, Bayesian)...

% Define a basis over the track...
basis.n = 30;   % number of basis functions
basis.s = 70;   % width parameter
[basis,yrbf] = get1Dbasis('vonmises',basis.n,rpos*pi,basis.s);

% Positions to use for decoding...
pvec = linspace(0,2*pi,256); % track length is mapped 0-2*pi for von-mises/fourier bases
[tmp,dbasis] = get1Dbasis('vonmises',basis.n,pvec,basis.s);

X = zscore(spk');
W = X\yrbf;
f = X*(W*dbasis');

% Transform and Standardize the read-out function...
f = exp(f);
f = bsxfun(@rdivide,f,sum(f,2));
decodePostProcess_basic
decoding_err(2,1) = mean(abs(err));
decoding_err(2,2) = median(abs(err));