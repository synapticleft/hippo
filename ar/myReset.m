function [x, V] = myReset(y, A, Q, R, init_x, init_V, spikes)
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
     [x(:,t), V(:,:,t)] = myUpdate(A, Q, R, y(:,t), prevx, prevV,initial,spikes(t));  
end

function [xnew, Vnew] = myUpdate(A, Q, R, y, x, V, initial,spike)
if initial
  xpred = x;
  Vpred = V;
else
  xpred = A*x;
  Vpred = A*V*A' + Q;
end
if spike
    xnew = y;
    Vnew = R;
else
   xnew = xpred;
  Vnew = Vpred;
end