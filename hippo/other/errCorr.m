function errCorr(v)

a = -arburg(v(:,1),1);
vHat = v(1:end-1,1)*a(2);
er = v(2:end,1) - vHat;
figure;plot(abs(xcov((er),(v(2:end,2)),1000,'coeff')));
%figure;plot(er);
figure;plot(abs(er));hold all;plot(abs(v(:,2))/10);
%[F0 S0]= arburg(v(:,1),1);
%F0 = F0(2);
%A0 = 1;
%W0 = 0.001;
%[A, W, F, S, ~, ~, x0, V0, loglik, xsmooth] = kalmanMLE(v(:,1), A0, W0, F0, S0, 0, 0);
%plot(real(v(:,1)));hold all;plot(real(xsmooth));