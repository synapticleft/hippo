function [A,W,F,S,xsmooth] = kalDeMod(X,Xf,inds)

X = X(:,inds);
Xf = Xf(:,inds);

%[u s v] = svds(Xf,1);
[L,Ph,LM,W0]=ffa(Xf.',1);
v = conj(Xf.'*pinv(L).');
[F0 S0]= arburg(v,1);%v(:,1)
F0 = -(F0(2));
A0 = L;
[A, W, F, S, ~, ~, x0, V0, loglik, xsmooth] = kalmanMLE(Xf.', A0, W0, F0, S0, 0, 0);
[F F0]
[S S0]
figure;subplot(121);imagesc(abs(W));subplot(122);imagesc(abs(W0))
%%for real-valued data
% A0 = [real(L) -imag(L)];
% x0 = [0;0];
% V0 = [S0 0;0 S0];
% F0 = [real(F0) -imag(F0); imag(F0) real(F0)];
% S0 = S0*ones(2);
% W0 = real(W0);
% [A, W, F, S, ~, ~, x0, V0, loglik, xsmooth] = kalmanMLE(real(X)', A0,W0, F0, S0, x0, V0);
% %figure;plot(xsmooth);hold all;plot([real(v) imag(v)])
% %return
% figure;plot(A0);hold all;plot(A);
% [F0 F]
% [S0 S]
% figure;subplot(121);imagesc(W);subplot(122);imagesc(W0)