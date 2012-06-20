% Slow.m - slowness penalty on coeff amplitudes
%

function S = SlowAR(Z,AR)
%i1 = I(:,2:end); I = I(:,1:end-1);
%AR = i1(:).'/I(:).';
S=Z(:,2:end)-AR*Z(:,1:end-1);
S = S.*conj(S);
