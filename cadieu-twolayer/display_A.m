%  display_network.m -- displays the state of the network 
%  (complex bfs and output variances)
%
%  h=display_network(A,Z_var,h);
%
%  A = complex basis function matrix
%  Z_var = vector of coefficient variances
%

function [array] =display_A(m,Z_var,fig_num)

if ~exist('fig_num','var')
    fig_num=1;
end
%m.e = [];
if isstruct(m)
    A = double(m.A);
else
    A = m;
    m.e = [];
end
if isfield(m,'dewhitenMatrix')
    A = m.dewhitenMatrix*A;
    array = display_Ahelper(A(1:end,:),fig_num);
    %A = m.zerophaseMatrix*A;
    display_Ahelper(m.zerophaseMatrix*A,fig_num+1);%A(1:end,:)
    if nargin < 2
        for i = 1:inf
            display_Ahelper(A*exp(1i*i/10),fig_num);
            display_Ahelper(m.zerophaseMatrix*A*exp(1i*i/10),fig_num+1);drawnow;
        end
    end
else
    display_Ahelper(A,fig_num);
end
if exist('Z_var','var')  && ~isempty(Z_var) && (max(Z_var(:)) > 0)
  sfigure(fig_num+2);
  subplot(211)
  bar(double(Z_var)), axis([0 m.N+1 0 double(max(Z_var))])
  title('Z variance')
  subplot(212)
  normA=sqrt(sum(abs(m.A).^2));
  plot(m.E);%bar(double(normA)), axis([0 m.N+1 0 double(max(normA))])
  title('basis norm (L2)')
end

drawnow

function array = display_Ahelper(A,fig_num)

[L M]=size(A);
sz = [8 L/8];
if sz(2) ~= round(sz(2))
    sz = [L 1];
end
if L == 96
    sz = [16 6];
end
%sz=sqrt(L);

buf=1;
m = ceil(sqrt(M));
n = ceil(M/m);

array=-ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf))*(1+1j);
k=1;
for c=1:m
  for r=1:n
      if k > M
          break;
      end
    clim=max(abs(A(:,k)));
    array(buf+(r-1)*(sz(1)+buf)+[1:sz(1)],buf+(c-1)*(sz(2)+buf)+[1:sz(2)])=...
	reshape(A(1:end,k),sz(1),sz(2))/clim;
    k=k+1;
  end
end

% if exist('h','var') && ~isempty(h)
%   subplot(221)
%   set(h(1),'CData',real(array));
%   axis image off
%   colormap gray; freezeColors
%   subplot(222)
%   set(h(2),'CData',imag(array));
%   axis image off
%   colormap gray; freezeColors
%   subplot(223)
%   set(h(3),'CData',abs(array));
%   axis image off
%   colormap gray; freezeColors
%   subplot(224)
%   set(h(4),'CData',angle(array));
%   alpha(abs(array)/max(abs(array(:))));
%   axis image off
%   colormap hsv; freezeColors
% else
subp_space = 0.03;
sfigure(fig_num);
clf;
%title('First-layer complex basis functions (m.A)')
colormap gray
subp(2,2,1,subp_space);
h(1)=imagesc(real(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('real')
subp(2,2,2,subp_space);
h(2)=imagesc(imag(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('imag')
subp(2,2,3,subp_space);
h(3)=imagesc(abs(array),[0 max(abs(array(:)))]);
axis image off
title('abs')
colormap gray; freezeColors
subp(2,2,4,subp_space);
h(4)=imagesc(angle(array),[-pi pi]);
alpha(abs(array)/max(abs(array(:))));
axis image off
colormap hsv; freezeColors
title('angle')
% end
