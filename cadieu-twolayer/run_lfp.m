function m = run_lfp(Xf)%pos,v,Xf,thresh)

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

% pos(pos == -1) = nan;
% reject = 0;
% for i = 1:4
%     reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
% end
% pos(reject,:) = nan;
% if size(v,1) < size(pos,1)
%     pos = pos(1:size(v,1),:);
% end
% for i = 1:4
%     nanInds = find(~isnan(pos(:,i)));
%     pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
% end
% nanInds = isnan(pos(:,1)) | isnan(pos(:,3));
% pos = pos(~nanInds,:);v = v(~nanInds,:);Xf = Xf(:,~nanInds);
% vel = angVel(pos);
% vel = [0; vel(:,1)];
% pos = bsxfun(@minus,pos,mean(pos));
% [a,~,~] = svd(pos(:,1:2),'econ');pos = a;
% for i = 1:2    
%     pos(:,i) = pos(:,i) - min(pos(:,i));
%     pos(:,i) = pos(:,i)/(max(pos(:,i)));
%     pos(:,i) = min(pos(:,i),.9999);
% end
% offSet = 1;
% Xf = [bsxfun(@times,Xf,exp(1i*angle(v(:,1).')));...
%      [zeros(offSet,1); v(1+offSet:end,1).*conj(v(1:end-offSet,1))./abs(v(1:end-offSet,1))].'];
% %Xf = [real(Xf);imag(Xf)];
% vel = filtLow(vel,1250/32,1);
% vel = vel/max(vel);inds = vel > thresh;

%% Initialize parameters

reset(RandStream.getDefaultStream);
warning('off','MATLAB:divideByZero')
warning('off','MATLAB:nearlySingularMatrix')

m.patch_sz =  size(Xf,1); % num elecs size
m.N = 256;  % firstlayer basis functions

% specify priors
p.firstlayer.prior = 'cauchy_Z';%'laplace_Z';%'l1l2';%'slow_cauchy';% changed per jascha's suggestion %slow_

% specify outerloop learning method
p.firstlayer.basis_method = 'steepest_adapt';%'steepest';%

% specifiy inference methods
p.firstlayer.inference_method='minFunc_ind';%'steepest';%

% misc
p.use_gpu = 0;
p.renorm_length=1;%1;
%p.normalize_crop=0;
p.whiten_patches=1;
p.p_every = 0;
p.show_p = 1;
p.quiet = 0;
%Xf = Xf(:,inds);
%% Init
[m, p] = init(m,p,Xf);%(:,inds)

% save path
fname=[strrep(datestr(now),' ','_') sprintf('patchsz%d_A%dx%d',m.patch_sz,m.M) '_%s.mat'];

% display parameters
display_every= 10;

%% learn firstlayer A
num_trials = 500;
save_every= num_trials;%100;
for i = 1:4
    learn_firstlayer
    p.firstlayer.eta_dA_target = p.firstlayer.eta_dA_target -.01;%/2;
end
