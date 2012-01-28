%
% function net=lds(X,K,cyc,tol);
% 
% Adaptive Linear Dynamical System 
%
% X - N x p data matrix
% K - size of state space (default 2)
% T - length of each sequence (N must evenly divide by T, default T=N)
% cyc - maximum number of cycles of EM (default 100)
% tol - termination tolerance (prop change in likelihood) (default 0.0001)
%
% net is a structure consisting of:
%
% A - state transition matrix
% C - observation (output) matrix 
% Q - state noise covariance 
% R - observation noise covariance
% x0 - initial state mean
% P0 - initial state covariance
% Mu - output mean
% LL - log likelihood curve
%
% Iterates until a proportional change < tol in the log likelihood 
% or cyc steps of EM

function [net net0]=gldsp(X,K,pOrder)

xform = 0;
rnd = 0;%.1;
if size(X,1) < size(X,2) X = X.'; end
[T,p] = size(X);

Mu=mean(X);
X=bsxfun(@minus,X,Mu);

    % initialize with Factor Analysis
    
    fprintf('\nInitializing with Factor Analysis...\n');
    [L,Ph,LM,Ph1]=ffa(X,K,100,0.001);
%    [u,~,~] = svds(X,K);
%    L = L*pinv(L)*u;
    R=Ph1;%diag(Ph);
    [A,~,Q] = mvar(X*pinv(L).',pOrder);
    [A,C,Q] = AR_to_SS(reshape(A,K,K,pOrder),Q(:,(end-K+1):end));
    C = L*C;
    x0 = zeros(K*pOrder,1);
    P0 = Q;
    %Phi=diag(1./R);
    fprintf('FA log likelihood %f\nInitialized.\n',LM(end));
    if xform
        [V,D] = eig(Q);
        Q = eye(size(Q,1));
        C = C*(V*sqrt(D));
        A = (sqrt(D)\V')*A*(V*sqrt(D));
    elseif rnd
        A = A.*(1+rnd*randn(size(A)));R = R.*(1+rnd*randn(size(R)));C = C.*(1+rnd*randn(size(C)));Q = Q.*(1+rnd*randn(size(Q)));
    end
net0.A = A;net0.R = R;net0.C = C;net0.Q = Q;net0.x0 = x0;net0.P0 = P0;
[C, R, A, Q, ~,~, x0, V0, loglik,xsmooth] = kalmanMLE(X,C,R,A,Q,x0,P0,0,0);
if xform
[V,D] = eig(Q);
Q = eye(size(Q,1));
C = C*(V*sqrt(D));
A = (sqrt(D)\V')*A*(V*sqrt(D));
end
net.C = C;net.R = R; net.A = A; net.Q = Q; net.x0 = x0; net.P0 = V0; net.loglik = loglik;net.xsmooth = xsmooth;