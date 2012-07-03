function A  = init_complex(M,N)
A=complex(randn(M,N),randn(M,N));
realnormA = sqrt(sum((real(A).^2)))';
imagnormA = sqrt(sum((imag(A).^2)))';
%if p.renorm_length
    % Keep our BFs from expanding; adapt magnitude of real and imaginary parts
    %realnormA = sqrt(sum((real(m.A).^2)))';
    %imagnormA = sqrt(sum((imag(m.A).^2)))';
    %m.A=complex(real(m.A)*diag(1./realnormA'),imag(m.A)*diag(1./imagnormA'));
    A = A*diag(1./sqrt(sum(A.*conj(A))));
%end
%A=complex(real(A)*diag(1./realnormA'),imag(A)*diag(1./imagnormA'));