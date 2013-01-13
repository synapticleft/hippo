function C = fastConv(A,B)

% perform "matrixized" fast convolution of matrix A and B along
% dimension dim. The convolution is fast because it uses the matrix
% operation (and fft/ifft) only. 
%
% for example:
%      fastConv(A,B,1) performs convolution of corresponding
%      columns in A,B
%      fastConv(A,B,2) performs convolution of corresponding
%      rows in A,B
%
% currently dim = 1 or 2
%
% the code is equivalent to running conv(A_i,B_i, 'full') in matlab
% (where A_i and B_i are columns (dim=1) or rows (dim=2) of A,B)
% and then stack the results together
%
% Zhen James Xiang
% Nov 14, Sun, 2010

% if 1==dim || nargin<3 % default
%   A = [A;zeros(size(A))];
%   B = [B;zeros(size(B))];
% elseif 2==dim
z1x = size(A,2);z2x = size(B,2);x1 = size(B,1);
% end
dim = 2;
C = ifft(fft([A,zeros(x1,z2x)],[],dim).*fft([B,zeros(x1,z1x)],[],dim),[],dim);

%z1x = size(A,2);z2x = size(B,2);
%C=real(ifft(fft(A,z1x+z2x-1,dim).*fft(B,z1x+z2x-1,dim),[],dim));%z1x+z2x-1
C = C(:,z2x:z1x);

%if 1==dim || nargin<3 % default
%  C = C(1:end-1,:);
%elseif 2==dim
%  C = C(:,1:end-1);
%end

