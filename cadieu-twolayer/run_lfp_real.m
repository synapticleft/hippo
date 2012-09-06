function run_lfp_real(pos,v,Xf,thresh)

%% driver script for learning a model
%%%   
%%%   m - struct containing the basis function variables
%%%     m.patch_sz - image domain patch size (side length of square patch)
%%%     m.M - whitened domain dimensions
%%%     m.N - firstlayer basis function dimensions
%%%     m.L - phase transformation basis functions dimensions
%%%     m.K - amplitude basis function dimensions
%%%     m.t - number of learning iterations
%%%     
%%%     m.A - first layer complex basis functions (m.M x m.N)
%%%     m.D - second layer transformation components (m.N x m.L)
%%%     m.B - second layer amplitude components (m.N x m.K)
%%%     
%%%   p - struct containing the learning, inference, and other parameters
%%%     p.firstlayer - first layer complex basis function parameters
%%%     p.ampmodel - second layer amplitude component parameters
%%%     p.phasetrans - second layer phase transformation parameters
%%%     
%%%   Summary of learning proceedure:
%%%       1. Initialize parameters
%%%       2. Estimate whitening transform
%%%       3. Learn first layer complex basis functions
%%%       4. Infer large batch of first layer coefficients
%%%       5. Learn second layer phase tranformation components
%%%       6. Learn second layer amplitude components
%%%       7. Display the results

bounds = [.1 .9];
pos(pos == -1) = nan;
reject = 0;
for i = 1:4
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
if size(v,1) < size(pos,1)
    pos = pos(1:size(v,1),:);
end
for i = 1:4
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end
nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);
vel = angVel(pos);
vel = [0; vel(:,1)];
pos = bsxfun(@minus,pos,mean(pos));
[a,~,~] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - min(pos(:,i));
    pos(:,i) = pos(:,i)/(max(pos(:,i)));
    pos(:,i) = min(pos(:,i),.9999);
end
offSet = 1;
Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1).')));...
     [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
Xf = [real(Xf);imag(Xf)];
Xf = bsxfun(@minus,Xf,mean(Xf,2));
vel = filtLow(vel,1250/32,1);
vel = vel/max(vel);inds = vel > thresh;

%% Initialize parameters

%% Whitening transform
p.whitening.pixel_noise_fractional_variance = .0002;
% the following seems redundant
% p.whitening.pixel_noise_variance_cutoff_ratio = 1; %1.25; % 1 + var(signal)/var(noise)
% p.whitening.X_noise_fraction = 8.;
% p.whitening.X_noise_var = .01;
C = Xf*Xf'/size(Xf,2);
[E, D] = eig(C);
[~, sind] = sort(diag(D),'descend');
d = diag(D);
d = d(sind);
E = E(:,sind);
% determine cutoff:
%variance_cutoff = p.whitening.pixel_noise_variance_cutoff_ratio*m.pixel_noise_variance;
rdim = size(Xf,1); % number of valid dims

factors = real(d(1:rdim).^(-.5));


noise_factors = ones(rdim,1);

%%removed per Jascha's suggestion

rolloff_ind = rdim-5;%rdim - sum(d(1:rdim)>variance_cutoff*p.whitening.X_noise_fraction);
noise_factors(rolloff_ind+1:end) = .5*(1+cos(linspace(0,pi-.01,rdim-rolloff_ind))); 
factors = factors.*noise_factors;
E = E(:,1:rdim);
D = diag(factors);
%% whitening transform
whitenMatrix = D*E';
dewhitenMatrix = E*D^(-1);
zerophaseMatrix = E*D*E';

Phi = sparsenet(whitenMatrix*Xf(:,inds),dewhitenMatrix,zerophaseMatrix);