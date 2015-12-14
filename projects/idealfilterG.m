function g = idealfilterG(f,pix,flip)

pFloor = floor(pix);
pCeil = ceil(pix);
pMod = mod(pix,1);
[M,N,P]=size(f);
u=0:(M-1);
v=0:(N-1);
idx=find(u>M/2);
u(idx)=u(idx)-M;
idy=find(v>N/2);
v(idy)=v(idy)-N;
[V,U]=meshgrid(v,u);
D=sqrt(U.^2+V.^2);
if exist('flip','var')
    H=(1-pMod)*double(D>pFloor)+pMod*double(D>pCeil);
else
    H=(1-pMod)*double(D<=pFloor)+pMod*double(D<=pCeil);
end
g = zeros(size(f));
for i = 1:P
    F = fft2(double(f(:,:,i)));
    G = H.*F;
    g(:,:,i) = real(ifft2(double(G)));
end