function array = showCompBases(A,sz)

[L M]=size(A);
if exist('sz','var')
    sz = [sz L/sz];
else
    sz = [8 L/8];
end

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
	reshape(A(:,k),sz(1),sz(2))/clim;
    k=k+1;
  end
end

imagesc(angle(array),[-pi pi]);
alpha(abs(array)/max(abs(array(:))));
axis image off
colormap hsv; freezeColors