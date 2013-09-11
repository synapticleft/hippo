function im = complexIm(in,sub,pow,scale,setMax,h)
in = squeeze(in);
if ~exist('scale','var') || isempty(scale)
    scale = 1;
end
if ~exist('sub','var') || isempty(sub)
    sub = 0;
end
if ~exist('pow','var') || isempty(pow)
    pow = 1;
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
    setMax = prctile(abs(in(:)),99.9);
end
if ~exist('h','var') || isempty(h)
    %im(:,:,1) = min(1,max(0,(angle(in)*scale)/(2*pi)+.5));
    %if scale > 1
    %    im(:,:,1) = (tanh(scale*angle(in)/pi)+1)/3;%mod(angle(in)*scale,2*pi)/(2*pi);
    %else
        im(:,:,1) = mod(angle(in)*scale,2*pi)/(2*pi);
   % end
    im(:,:,2) = 1;im(:,:,3) = min(1,power(abs(in)/setMax,pow));
else
    im(:,:,1) = angle(h);im(:,:,3) = min(1,abs(h)./setMax);im(:,:,2) = min(1,power(abs(in)/setMax,pow));
end
temp = squeeze(im(:,:,3));
temp(isnan(in)) = 1;
im(:,:,3) = temp;
temp = squeeze(im(:,:,1));
temp(isnan(in)) = 1;
im(:,:,1) = temp;
temp = squeeze(im(:,:,2));
temp(isnan(in)) = 0;
im(:,:,2) = temp;
%im(:,:,3) = max(squeeze(im(:,:,3)),isNan);
im = hsv2rgb(im);