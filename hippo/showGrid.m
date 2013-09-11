function array=showGrid(A,grid,rad,ratio,in,myDim,clim)
%% image in a grid data sets that have multiple images
%% like electrode arrays or activations of various ICs

if ~exist('ratio','var') || isempty(ratio)
    ratio = [1 1];
end

buf=1;
%if ndims(A) == 2 %&& numel(A) ~= prod(size(A))
%    A = A.';
%end
if exist('myDim','var') && ~isempty(myDim)
    n = myDim; m = ceil(size(A,1)/n);
else
n = ceil(sqrt(size(A,1)));
m = ceil(size(A,1)/n);
end
if ndims(A) == 2
if numel(grid) == 2
    sz = grid;
    grid = reshape(1:prod(grid),sz);
else
    sz = size(grid);
    grid = grid-min(grid(:)) + 1;
end
max(grid(:))
A1 = zeros(size(A,1),sz(1),sz(2));
A1(:,:) = A(:,grid(:));
else
    A1 = A;
    sz = size(A1);sz = sz(2:3);
end
if exist('in','var') && ~isempty(in) && in
    [xout yout] = meshgrid(ratio(2):ratio(2)*sz(2),ratio(1):ratio(1)*sz(1));
    [xin yin] = meshgrid((1:sz(2))*ratio(2),(1:sz(1))*ratio(1));
    sz = size(xout);
end

array=ones(buf+m*(sz(1)+buf),buf+n*(sz(2)+buf));
if ~isreal(A)
    array = array*nan;
end

%k=1;
%quivDat = zeros(size(A,2),4);
for j=1:m
  for i=1:n
      k = (j-1)*n+i;
    if k > size(A,1)
        break
    end
    %temp = reshape(A(:,k),sz(1),sz(2))/clima;
    indsx = buf+(j-1)*(sz(1)+buf)+[1:sz(1)];
    indsy = buf+(i-1)*(sz(2)+buf)+[1:sz(2)];
%    [fx fy] = gradient(temp);
    if exist('in','var') && ~isempty(in) && in
        temp = interp2(xin,yin,squeeze(A1(k,:,:)),xout,yout,'cubic');
        %array(indsx,indsy) = temp;
    else
        temp = squeeze(A1(k,:,:));%temp;
    end
    if exist('rad','var') && ~isempty(rad) && rad
        temp = imfilter(temp,fspecial('gaussian',max(5,rad*2),rad),'replicate');
    end
    if isreal(temp)
        temp = temp - prctile(temp(:),5);%min(temp(:));
    end
    if ~exist('clim','var') 
        clima=max(abs(temp(:)));
    else
        clima = clim;
    end
    %temp = temp/max(temp(:));
   % temp = rot90(rot90((temp)));
    array(indsx,indsy) = temp/clima;
%    quivDat(k,:) = [mean(indsx) mean(indsy) mean(fx(:)) mean(fy(:))];
%    k=k+1;
  end
end
%array = imfilter(array,fspecial('gaussian',5,.01));
%scale = 10;
%colormap gray
%subplot(211)
array = array(2:end-1,2:end-1);
x = 1:size(array,2);
y = 1:size(array,1);
if ~exist('in','var') || (~isempty(in) && ~in)
x = x*ratio(1);y = y*ratio(2);
end
%if ~exist('h','var') || isempty(h)
%    if nargout>0
%        hout=imagesc(x,y,array,[-1 1]);
%    else
%array = imfilter(array,fspecial('gaussian',5,1));

if isreal(array)
    imagesc(x,y,array,[0 1]);%[min(array(array~=-1)) max(array(array~=-1))]);colormap gray;%
    colormap gray;
else
    imagesc(x,y,complexIm(array,0,.5,[],1));
end
if exist('in','var') && ~isempty(in) && in
    axis image;
end
%red = 6;
%array1 = imresize(array,1/red,'bilinear');
% [dx dy] = angGradient(array1);
% [xs ys] = meshgrid(red*(1:size(array1,2))-2,red*(1:size(array1,1))-2);
% hold all;quiver(xs,ys,-dx,-dy,'color','w','linewidth',2);
axis off;
%axis image
%    end
%    hold on;
%    quiver(quivDat(:,2),quivDat(:,1),quivDat(:,3),quivDat(:,4),'r','linewidth',3);
%    hold off;
%    axis off
%else
%    set(h,'CData',array)
%end
% s
% subplot(212)
% 
% normA=sqrt(sum(A.*A));
% bar(normA), axis([0 M+1 0 max(normA)])
% title('basis norm (L2)')


    drawnow