% Make a point move in the 2D plane
% State = (x y xdot ydot). We only observe (x y).
% Generate data from this process, and try to learn the dynamics back.
% HIGHER ORDER
% X(t+1) = A X(t) + noise(Q)
% Y(t) = C X(t) + noise(R)

ss = 2;p = 2;
[u s v] = svds(Xf(:,1:1000),ss);
[A,Q] = mvar(v,p);
[A,C,Q,R] = AR_to_SS(reshape(A,ss,ss,p),Q(:,(end-ss+1):end));
R = rand(ss).*eye(ss);
initx = .1*ones(1,ss);
T = 1000;
%[x,y] = sample_lds(A, C, Q, R, initx, T);
figure;subplot(211);plot(real(x)');
[net,net0] = gldsp(Xf(:,1:1000),ss,p);
subplot(212);plot(real(net.xsmooth));
% Initializing the params to sensible values is crucial.
% Here, we use the true values for everything except F and H,
% which we initialize randomly (bad idea!)
% Lack of identifiability means the learned params. are often far from the true ones.
% All that EM guarantees is that the likelihood will increase.
%F1 = randn(ss,ss);
%H1 = randn(os,ss);
%Q1 = Q;
%R1 = R;
%initx1 = initx;
%initV1 = initV;
%max_iter = 10;
%[F2, H2, Q2, R2, initx2, initV2, LL] =  learn_kalman(y, F1, H1, Q1, R1, initx1, initV1, max_iter);

