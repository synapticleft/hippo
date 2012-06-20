% Slowp.m - slowness penalty derivative
%

function Sp = SlowARp(Z,AR)
%i1 = I(:,2:end); I = I(:,1:end-1);
%AR = i1(:).'/I(:).';
D=Z(:,2:end)-AR*Z(:,1:end-1);
%D=diff(a,1,2);

Sp=[-D(:,1)*conj(AR) -(D(:,2:end)*conj(AR)-D(:,1:end-1)) D(:,end)];
