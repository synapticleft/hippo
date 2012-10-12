function im = complexIm(in,sub,pow,scale,h,setMax)
if ~exist('scale','var') || isempty(scale)
    scale = 1;
end
if sub
    s = size(in);
    in = in(:);
    in = in-mean(in);
    in = [real(in) imag(in)];
    in = in/sqrtm(cov(in));
    in  = reshape(complex(in(:,1),in(:,2)),s);
    %in = in-mean(in(:));
end

if ~exist('setMax','var')
    setMax = prctile(abs(in(:)),99);
end
if ~exist('h','var') || isempty(h)
    %im(:,:,1) = min(1,max(0,(angle(in)*scale)/(2*pi)+.5));
    im(:,:,1) = mod(angle(in)*scale,2*pi)/(2*pi);
    im(:,:,2) = 1;im(:,:,3) = min(1,power(abs(in)/setMax,pow));
else
    im(:,:,1) = angle(h);im(:,:,3) = min(1,abs(h)./setMax);im(:,:,2) = min(1,power(abs(in)/setMax,pow));
end
im = hsv2rgb(im);