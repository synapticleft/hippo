function [Ierror, Ihat] = calc_Ierror(I,Z,m,p)%, ampError, phaseError

%Z = a.*exp(1j*phase);
%Ihat = real(m.A)*(a.*cos(phase)) + imag(m.A)*(a.*sin(phase));%real(Phi(:,:,t)*conj(Z));
Ihat = m.A*conj(Z);%
%figure(8);plot(sum(a));hold all;plot(real(sum(Ihat)),'r');hold off;drawnow;
Ierror = I - Ihat;
% ampError = abs(I) - abs(Ihat);
% phaseError = angle(I) - angle(Ihat);
% phaseError = mod(phaseError+pi,2*pi)-pi;