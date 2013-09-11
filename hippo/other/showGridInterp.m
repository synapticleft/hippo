function hout = showGridInterp(A,original,sz,clim)
%%artifact of showGrid

[L M]=size(A);
buf=1;
    n = ceil(sqrt(size(A,2)));
    m = ceil(size(A,2)/n);
    if numel(sz) == 1
        sz = [sz sz];
    end
    array=-ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf));

k=1;
mn = min(original);
mx = max(original);
for j=1:m
  for i=1:n
    if k > size(A,2)
        break
    end
    if ~exist('clim','var') 
        clima=max(abs(A(:,k)));
    else
        clima = clim;
    end
    [size(original) size(A)]
size(interp2(original(1,:),original(2,:),A(:,k),linspace(mn(1),mx(1),sz(1)),linspace(mn(2),mx(2),sz(2)))/clima);
return
    array(buf+(i-1)*(sz(1)+buf)+[1:sz(1)],buf+(j-1)*(sz(2)+buf)+[1:sz(2)])=...
	reshape(A(:,k),sz(1),sz(2))/clima;
    k=k+1;
  end
end

colormap gray
%subplot(211)
x = 1:size(array,1);
y = 1:size(array,2);
x = x*ratio(1);y = y*ratio(2);
if ~exist('h','var') || isempty(h)
    if nargout>0
        hout=imagesc(x,y,array,[-1 1]);
    else
        imagesc(x,y,array,[-1 1])
    end
    axis off
else
    set(h,'CData',array)
    drawnow
end