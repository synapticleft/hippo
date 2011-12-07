% showbfs.m - function to show basis functions
%
% function hout=showGrid(A,bg,h)
%
% A = bf matrix
% bg = 'black' or 'white' (default='black')

function hout=showJru(A,gridsize,ratio)

bg='black';
clim = 1;
[L M]=size(A);
buf=1;

%if size(A,2) == prod(gridsize)
%    n = gridsize(1);m = gridsize(2);
%else
    n = size(A,1);%ceil(sqrt(size(A,2)));
    m = size(A,2);%ceil(size(A,2)/n);
%end
sz = gridsize;

if bg=='black'
  array=-ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf));
else
  array=ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf));
end

for j=1:m
  for i=1:n
%     if ~exist('clim','var') 
%         clim=max(abs(A(:,k)));
%     end
    array(buf+(i-1)*(sz(1)+buf)+[1:sz(1)],buf+(j-1)*(sz(2)+buf)+[1:sz(2)])=...
	reshape(A(j,i,:),sz(1),sz(2))/clim;
  end
end

x = 1:size(array,1);
y = 1:size(array,2);
x = x*ratio(1);y = y*ratio(2);
if nargout>0
    hout=imagesc(x,y,array,[-1 1]);
else
    imagesc(x,y,array,[-1 1])
end
axis image off
drawnow

