function im = complexIm(in,sub,pow,scale,h)
if ~exist('scale','var')
    scale = 1;
end
if sub
    in = in-mean(in(:));
end
if ~exist('h','var')
    im(:,:,1) = min(1,max(0,(angle(in)*scale)/(2*pi)+.5));im(:,:,2) = 1;im(:,:,3) = min(1,power(abs(in)/prctile(abs(in(:)),99),pow));
else
    im(:,:,1) = angle(h);im(:,:,3) = min(1,abs(h)./prctile(abs(h(:)),99));im(:,:,2) = min(1,power(abs(in)/prctile(abs(in(:)),99),pow));
end
im = hsv2rgb(im);