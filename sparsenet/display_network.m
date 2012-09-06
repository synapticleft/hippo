function display_network(A,S_var)
%
%  display_network -- displays the state of the network (weights and 
%                     output variances)
%
%  Usage:
%
%    display_network(A,S_var);
%
%    A = basis function matrix
%    S_var = vector of coefficient variances
%

figure(1)

[L M]=size(A);

%sz=sqrt(L);
sz = [8 4];
buf=1;

if floor(sqrt(M))^2 ~= M
  n=sqrt(M/2);
  m=M/n;
else
  m=sqrt(M);
  n=m;
end

array=-ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf));

k=1;

for j=1:m
  for i=1:n
    clim=max(abs(A(:,k)));
    array(buf+(i-1)*(sz(1)+buf)+[1:sz(1)],buf+(j-1)*(sz(2)+buf)+[1:sz(2)])=...
	reshape(complex(A(1:32,k),A(34:65,k)),sz(1),sz(2))/clim;
    k=k+1;
  end
end

subplot(212)
imagesc(array,[-1 1]), axis image off
title('basis functions')

if exist('S_var','var')
  subplot(211)
  bar(S_var), axis([0 M+1 0 max(S_var)])
  title('coeff variance')
end

drawnow
