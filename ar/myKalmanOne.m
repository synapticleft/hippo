function [x, V, VV,loglik] = myKalmanOne(y, A, Q, R, init_x, init_V, spikes)
% INPUTS:
% y(:,t)   - the observation at time t
% A - the system matrix
% Q - the system covariance 
% R - the observation covariance
% init_x - the initial state (column) vector 
% init_V - the initial state covariance 
% OUTPUTS (where X is the hidden state being estimated)
% x(:,t) = E[X(:,t) | y(:,1:t)]
% V(:,:,t) = Cov[X(:,t) | y(:,1:t)

[~, T] = size(y);
ss = size(A,1); % size of state space
x = zeros(ss, T);
V = zeros(ss, ss, T);
VV= zeros(ss, ss, T);
R = R(1,1);
C = [1 zeros(1,ss-1)];%[1 0 0 0];
y = y(1,:);
loglik = 0;
for t=1:T
  if t==1
    prevx = init_x;
    prevV = init_V;
    initial = 1;
  else
    prevx = x(:,t-1);
    prevV = V(:,:,t-1);
    initial = 0;
  end
    [x(:,t), V(:,:,t), VV(:,:,t), LL] = ..., 
	myUpdate(A, C, Q, R, y(t), prevx, prevV,initial,spikes(t));%,init_V
  loglik = loglik + LL;
end

function [xnew, Vnew, VVnew,loglik] = myUpdate(A, C, Q, R, y, x, V, initial,spike)
%y = y*(.5+randn/1000);
if initial
  xpred = x;
  Vpred = V;
else
  xpred = A*x;% + sqrt(Q(:,1))*complex(randn,randn)/sqrt(2);
  Vpred = A*V*A' + Q;
%  Vpred(1,1) = real(Vpred(1,1));
end
e = y - C*xpred; % error (innovation)
ss = length(V);
S = C*Vpred*C' + R;
loglik = 0;%myGprob(abs(e), zeros(1,length(e)), S);
if spike
    %Sinv = inv(S);
    K = Vpred*C'/S;%*Sinv; % Kalman gain matrix
else
    %S = Vpred + zSig*100000;
    %K = Vpred/S;
    K = zeros(ss,1);
end
% If there is no observation vector, set K = zeros(ss).
% if spike
%     xnew = y;
%     Vnew = R;
% else
   xnew = xpred + K*e;
%   Vnew = Vpred;
% end
%Vnew = (eye(ss) - K*C)*Vpred; %Old way, has instability
Vnew = (eye(ss)-K*C)*Vpred*(eye(ss)-K*C)'+K*R*K';
VVnew = (eye(ss) - K*C)*A*V;
