function [Z Ierror exit_flag]= infer_Z(I,m,p)
% infer_Z.m - infer firstlayer latent variables

sz = size(I,2);
% Initialize the latent variables
if p.use_gpu
    Z = .2*complex(grandn(m.N,sz),grandn(m.N,sz));
else
    Z = .2*complex(randn(m.N,sz),randn(m.N,sz));
%      Z = pinv(m.A)*I;
%      Z = conj(Z);
end
a = abs(Z);
phase = angle(Z);
astop = sz*m.N;
% aphase0 = [reshape(a,numel(a),1); reshape(phase,numel(phase),1)];
Z0 = [reshape(real(Z),numel(Z),1); reshape(imag(Z),numel(Z),1)];%reshape(Z,numel(Z),1);

% Setup parameters for the specified method
switch p.firstlayer.inference_method
    case 'steepest'
        % Generative image model
        %Ih = real(A*conj(Z));
        [E0, ~, ~, ~] = obj_fun_z(aphase0,I,m,p);
        exit_flag=1;
        for t=1:p.firstlayer.iter
            aphase0 = [reshape(a,numel(a),1); reshape(phase,numel(phase),1)];
            [E, daphase, Ih, Ierror] = obj_fun_z(aphase0,I,m,p);
            da = reshape(daphase(1:astop),m.N,sz);
            dphase = reshape(daphase((astop+1):2*astop),m.N,sz);
            dE=(E0-E)/E0;
            if (dE< -.01) && (t>10)
                exit_flag=0;
                fprintf('\rInference unstable... exiting')
                break
            elseif (dE<.00001) && (t>(.5*p.firstlayer.iter))
                exit_flag=1;
                fprintf('\rConverged at inter #: %1i',t)
                break
            end
            E0=E;
            % Update a, phase
            a     = a-p.firstlayer.eta_a*da;
            phase = phase-p.firstlayer.eta_phase*dphase;
            % deal with negative a
            anegind=a<=0;
            a(anegind)=1;
%            display_infer_Z(a,phase,I,Ih,m,p)
%             %phase(anegind)=angle(Zr(anegind));
%             Zr_real = real(m.A).'*Ierror;
%             Zr_imag = imag(m.A).'*Ierror;
%             Zr = complex(Zr_real,Zr_imag);
%             aneg_angle = angle(Zr);
%             phase(anegind) = 0;
%             phase = phase + anegind.*aneg_angle;
        end
        fprintf('.\n')
    case {'minFunc_ind','minFunc_ind_lbfgs'}
        [E, ~, ~, Ierror] = obj_fun_z(Z0,I,m,p);
        SNR = -10*log10(var(Ierror(:))/var(I(:)));        
        fprintf('\rE=%02.4e, SNR=%2.2f',double(E),double(SNR));
        [ZF, E, ~] = minFunc_ind(@obj_fun_z,Z0,p.firstlayer.minFunc_ind_Opts,I,m,p);
        Z = reshape(complex(ZF(1:astop),ZF((astop+1):end)),m.N,sz);
%         a = reshape(aphase(1:astop),m.N,sz);
%         phase = reshape(aphase((astop+1):2*astop),m.N,sz);      
        [Ierror, Ih] = calc_Ierror(I,Z,m,p);
        SNR = -10*log10(var(Ierror(:))/var(I(:)));
        fprintf('\rE=%02.4e, SNR=%2.2f\r\n',double(E),double(SNR));
        exit_flag=1;
end
if p.show_p
    display_infer_Z(Z,I,Ih,m,p)
end
%Z = a.*exp(1j*phase);
