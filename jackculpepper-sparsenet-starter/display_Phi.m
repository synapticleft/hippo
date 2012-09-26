function array = display_Phi(A,fig_num)

[L M]=size(A);
sz = [8 L/8];
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
	reshape(A(:,k),sz(1),sz(2))/clim;
    k=k+1;
  end
end

subp_space = 0.03;
sfigure(fig_num);
clf;
colormap gray
subp(2,3,1,subp_space);
h(1)=imagesc(real(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('real')
subp(2,3,2,subp_space);
h(2)=imagesc(imag(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('imag')
subp(2,3,4,subp_space);
h(3)=imagesc(abs(array),[0 max(abs(array(:)))]);
axis image off
title('abs')
colormap gray; freezeColors
subp(2,3,5,subp_space);
h(4)=imagesc(angle(array),[-pi pi]);
alpha(abs(array)/max(abs(array(:))));
axis image off
colormap hsv; freezeColors
title('angle')
drawnow;