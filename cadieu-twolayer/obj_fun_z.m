function [E, g, Ihat, Ierror] = obj_fun_z(X,I,m,p)

sz = size(I,2);
astop = sz*m.N;
Z = reshape(complex(X(1:astop),X((astop+1):end)),m.N,sz);
%Zimag = reshape(,m.N,sz);
%Z = complex(Zreal,Zimag);
%Z = reshape(X,m.N,sz);%(1:astop)

[Ierror, Ihat] = calc_Ierror(I,Z,m,p);
% Compute Energy Terms
switch p.firstlayer.prior
    case 'cauchy'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_cauchy(a(:),p.firstlayer.a_cauchy_beta,p.firstlayer.a_cauchy_sigma));
        E= mse + a_sparsity;
    case 'gauss'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_gauss(a(:),p.firstlayer.a_gauss_beta));
        E= mse + a_sparsity;
    case 'laplace'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_laplace(a(:),p.firstlayer.a_laplace_beta));
        E= mse + a_sparsity;
    case 'jascha'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_jascha(a(:),p));
        a_slowness = .5*p.firstlayer.a_lambda_S*sum(sum((Slow(a))));
        E= mse + a_sparsity + a_slowness;
    case 'slow_cauchy'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_cauchy(a(:),p.firstlayer.a_cauchy_beta,p.firstlayer.a_cauchy_sigma));
        a_slowness = .5*p.firstlayer.a_lambda_S*sum(sum((Slow(a))));
        E= mse + a_sparsity + a_slowness;
    case 'slow_gauss'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_gauss(a(:),p.firstlayer.a_gauss_beta));
        a_slowness = .5*p.firstlayer.a_lambda_S*sum(sum((Slow(a))));
        E= mse + a_sparsity + a_slowness;
    case 'slow_laplace'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_laplace(a(:),p.firstlayer.a_laplace_beta));
        slowCoeff = Slow(a);%gSlow(a,p);
        a_slowness = .5*p.firstlayer.a_lambda_S*sum(slowCoeff(:)); %
        E= mse + a_sparsity + a_slowness;
    case 'der_laplace'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_laplace(a(:),p.firstlayer.a_laplace_beta));
        sparseDer = p.firstlayer.a_lambda_S*SlowDer(a);
        slowTheta = .5*p.firstlayer.p_lambda_S*thetaSlow(phase);%
        E= mse + a_sparsity + sum(sparseDer(:)) + sum(slowTheta(:));
    case 'laplace_Z'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_laplaceZ(Z,p.firstlayer.a_laplace_beta));
        if p.firstlayer.a_lambda_S
            a_sparsity = a_sparsity + .5*p.firstlayer.a_lambda_S*sum(sum((Slow(abs(Z)))));
        end
        E= mse + a_sparsity;% + a_slowness;
    case 'laplace_AR'
        mse = sum(sum(bsxfun(@times,0.5*m.I_noise_factors,Ierror.*conj(Ierror))));
        a_sparsity = sum(S_laplaceZ(Z(:),p.firstlayer.a_laplace_beta));
        if p.firstlayer.a_lambda_S
            a_sparsity = a_sparsity + .5*p.firstlayer.a_lambda_S*sum(sum(SlowAR(Z,m.AR)));
        end
        E= mse + a_sparsity;% + a_slowness;
end
E = double(E);
if nargout>1
    weight_dEdz = conj(m.A'*bsxfun(@times,-m.I_noise_factors,Ierror));%
    %cosPhase = cos(phase);sinPhase = sin(phase);
    %grada = cosPhase.*real(weight_dEdz) + sinPhase.*imag(weight_dEdz);
    %gradphase = a.*(cosPhase.*imag(weight_dEdz) - sinPhase.*real(weight_dEdz));%
    switch p.firstlayer.prior
        case 'cauchy'
            grada = grada + dS_cauchy(a,p.firstlayer.a_cauchy_beta,p.firstlayer.a_cauchy_sigma);
        case 'gauss'
            grada = grada + dS_gauss(a,p.firstlayer.a_gauss_beta);
        case 'laplace'
            grada = grada + dS_laplace(a,p.firstlayer.a_laplace_beta);
        case 'jascha'
            grada = grada + dS_jascha(a,p)+ p.firstlayer.a_lambda_S*Slowp(a);
        case 'slow_cauchy'
            grada = grada + dS_cauchy(a,p.firstlayer.a_cauchy_beta,p.firstlayer.a_cauchy_sigma) ...
                          + p.firstlayer.a_lambda_S*Slowp(a);
        case 'slow_gauss'
            grada = grada + dS_gauss(a,p.firstlayer.a_gauss_beta) ...
                          + p.firstlayer.a_lambda_S*Slowp(a);
        case 'slow_laplace'
            grada = grada + dS_laplace(a,p.firstlayer.a_laplace_beta) ...
                          + p.firstlayer.a_lambda_S*Slowp(a);%(a-slowCoeff); %gSlowp(a,p) %
        case 'der_laplace'
            grada = grada + dS_laplace(a,p.firstlayer.a_laplace_beta) ...
                          + dSlowDer(a)*p.firstlayer.a_lambda_S;%
            gradphase = gradphase + dthetaSlow(phase)*p.firstlayer.p_lambda_S;
        case 'laplace_Z'
            absGrad = dS_laplaceZ(Z,p.firstlayer.a_laplace_beta);
            gradZ = weight_dEdz + [p.firstlayer.a_laplace_beta(1)*absGrad(1,:); ...
                p.firstlayer.a_laplace_beta(end)*absGrad(2:end,:)];%p.firstlayer.a_laplace_beta*
            if p.firstlayer.a_lambda_S
                gradZ = gradZ + p.firstlayer.a_lambda_S*Slowp(abs(Z)).*absGrad;
            end
        case 'laplace_AR'
            absGrad = dS_laplaceZ(Z,p.firstlayer.a_laplace_beta);            
            gradZ = weight_dEdz + [p.firstlayer.a_laplace_beta(1)*absGrad(1,:); ...
                p.firstlayer.a_laplace_beta(end)*absGrad(2:end,:)];%p.firstlayer.a_laplace_beta*
            if p.firstlayer.a_lambda_S
                gradZ = gradZ + p.firstlayer.a_lambda_S*SlowARp(Z,m.AR);%.*absGrad;
            end
    end
    g = [reshape(real(gradZ),numel(gradZ),1); reshape(imag(gradZ),numel(gradZ),1)];
end

%old attempts
%    weighted_error = bsxfun(@times,-m.I_noise_factors,Ierror);
%     %grada = (real(m.A).'*weighted_error).*cos(phase) + (imag(m.A).'*weighted_error).*sin(phase);
%     grada = real((m.A.'*weighted_error).*conj(exp(1i*phase)));
%     if p.firstlayer.natural_gradient
% %        gradphase = (imag(m.A).'*weighted_error).*cos(phase) - (real(m.A).'*weighted_error).*sin(phase);
%         gradphase = imag((m.A.'*weighted_error).*conj(exp(1i*phase)));
%     else
% %        gradphase = a.*((imag(m.A).'*weighted_error).*cos(phase) - (real(m.A).'*weighted_error).*sin(phase));
%         gradphase = a.*imag((m.A.'*weighted_error).*conj(exp(1i*phase)));
%     end

%dataProj = m.A'*I;
%modelProj = m.A'*m.A*(Z);%conj
% grada = abs(dataProj)-abs(modelProj);%abs(weighted_error);
% gradphase = angle(dataProj)-angle(modelProj);%angle(weighted_error);
% gradphase = mod(gradphase+pi,2*pi)-pi;
% grada = bsxfun(@times,-m.I_noise_factors,grada);
% gradphase = bsxfun(@times,-m.I_noise_factors,gradphase);%.*a