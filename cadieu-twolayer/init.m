function [m,p] = init(m,p,Xf)

%% misc
m.t = 0;
p.imszt=256;% number of time steps

%% whitening params %%

if p.whiten_patches
    p.whitening.pixel_noise_fractional_variance = .001;
    % the following seems redundant
    p.whitening.pixel_noise_variance_cutoff_ratio = 1; %1.25; % 1 + var(signal)/var(noise)
    p.whitening.X_noise_fraction = 5.;
    p.whitening.X_noise_var = .01;
    % run whitening
    p.whitening.whiten_num_patches = 20000;%min(400*m.patch_sz,200000)/20;%160000; TEMPORARY /20
    [m, p] = learn_whitening(m,p,Xf);
    %m.N = m.M;
else
    m.M = m.patch_sz;%m.N;
    m.I_noise_factors = 100*ones(m.M,1);%
    m.imageMean = mean(Xf,2);
    m.imageStd = std(Xf,0,2);
    %p.var = sqrt(10*var(Xf(:)));
end
%% init basis functions %%
m.A = init_complex(m.M,m.N);

%% first layer %%
p.firstlayer.use_GS = 1;
switch p.firstlayer.basis_method
    case 'steepest_adapt'
        p.firstlayer.A_eta=.01;
        p.firstlayer.eta_dA_target = .1;
        p.firstlayer.up_factor = 1.02;
        p.firstlayer.down_factor = .95;
        
    case 'steepest'
        p.firstlayer.A_eta=.05;
end

switch p.firstlayer.prior
    case 'slow_cauchy'        
        % a
        p.firstlayer.a_cauchy_beta = 10; % 2.2%1;%
        p.firstlayer.a_cauchy_sigma = .4; % .1
        p.firstlayer.a_lambda_S= .5;%5;
        p.firstlayer.a_thresh  = exp(-4);
    case 'cauchy'        
        % a
        p.firstlayer.a_cauchy_beta = 1;%10; % 2.2%1;%
        p.firstlayer.a_cauchy_sigma = .4; % .1
    case 'laplace'
        p.firstlayer.a_laplace_beta = 1000;
    case 'l1l2'
        p.firstlayer.a_laplace_beta = 2000;
    case 'slow_laplace'
        p.firstlayer.a_laplace_beta = [10 1000];%1;
        p.firstlayer.a_lambda_S = 10000;%.5;
        p.firstlayer.a_tau_S = 100;
    case 'laplace_Z'
        p.firstlayer.a_laplace_beta = 1000;
        p.firstlayer.a_lambda_S =0;
        %p.firstlayer.a_tau_S = 100;
    case 'laplace_AR'
        p.firstlayer.a_laplace_beta = 100*[0.99 1];%1;
        p.firstlayer.a_lambda_S = 0;%.5;
end

switch p.firstlayer.inference_method
    case 'steepest'
        p.firstlayer.iter  =  120;
        p.firstlayer.eta_a     = .005;%.00005
        p.firstlayer.eta_phase = .0005;
        p.firstlayer.natural_gradient = 1;
        
    case 'minFunc_ind'
        p.firstlayer.minFunc_ind_Opts.Method = 'bb';%'cg';%'bb';%'csd';%
        p.firstlayer.minFunc_ind_Opts.Display = 'off';
        p.firstlayer.minFunc_ind_Opts.MaxIter = 30;%15;%
        p.firstlayer.minFunc_ind_Opts.MaxFunEvals = 60;%20;%
        p.firstlayer.natural_gradient = 1;%1 originally

end
