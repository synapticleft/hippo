function [a,a2] = aarg(v,Xf,UC)
% Calculates adaptive autoregressive (AAR) and adaptive autoregressive moving average estimates (AARMA)
% of real-valued data series using Kalman filter algorithm.
% The AAR process is described as following  
%       y(k) - a(k,1)*y(t-1) = e(k);
% Input:
%       y       Signal (AR-Process)
%       Mode    aMode=1, LMS algorithm (signal normalized)
%               aMode=2, LMS algorithm with adaptive normalization              
%       MOP     model order, default [10,0] 
%               MOP=[p]         AAR(p) model. p AR parameters
%       UC      Update Coefficient, default 0     
% Output:
%       a       AAR(MA) estimates [a(k,1), a(k,2), ..., a(k,p),b(k,1),b(k,2), ..., b(k,q]
%       e       error process (Adaptively filtered process)
v = v/std(v);Xf = Xf/std(Xf(:));
[a] = getCoeffs([0; v(1:end-1)],v,UC);
[a2] = getCoeffs(v(1:end),Xf(:,1:end)',UC);

function a = getCoeffs(X,y,UC)
MOP = 1;aMode = 1;
[nc,~]=size(y);

a0=zeros(1,MOP);
MSY=mean(abs(y(:)).^2);

e=zeros(nc,1);
V=zeros(nc,1);
a=zeros(size(y));%a0(ones(nc,1),:);

%------------------------------------------------
%       First Iteration
%------------------------------------------------
Y=zeros(MOP,1);
E=y(1,:);
%e(1)=E;
V(1) = (1-UC) + UC*E*E';
%------------------------------------------------
%       Update Equations
%------------------------------------------------
for t=1:nc,
    if t == 1
        prevA = a0;
    else
        prevA = a(t-1,:);
    end
    % Prediction Error
    E = y(t,:) - prevA*X(t);
%    e(t) = E;
    E2=E'*E;
    if aMode == 1, % LMS
        %       V(t) = V(t-1)*(1-UC)+UC*E2;
        a(t,:)=prevA + ((UC/MSY)*E.'*X(t)').';
    elseif aMode == 2, % LMS with adaptive estimation of the variance
        V(t) = V(t-1)*(1-UC)+UC*E2;
        a(t,:)=prevA + UC/V(t)*E*X(t)';
    end
end
for t=nc-1:-1:1
    prevA = a(t+1,:);
    % Prediction Error
    E = y(t,:) - prevA*X(t);
%    e(t) = E;
    E2=E'*E;
    if aMode == 1, % LMS
        %       V(t) = V(t-1)*(1-UC)+UC*E2;
        a(t,:)=prevA + ((UC/MSY)*E.'*X(t)').';
    elseif aMode == 2, % LMS with adaptive estimation of the variance
        V(t) = V(t-1)*(1-UC)+UC*E2;
        a(t,:)=prevA + UC/V(t)*E*X(t)';
    end
end