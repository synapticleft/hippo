function sig = filtHigh(sig,Fs,f,order)

%used X = filtHigh(X,1250/8,4,4);

if nargin < 4
    order = 2;
end
nanInds = isnan(sig);
sig(nanInds) = 0;
f1 = fdesign.highpass('n,f3dB',order,f,Fs);
d1 = design(f1,'butter');
%a1 = filter(d1,flipud(filter(d1,flipud(a'))))';
[B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
if isreal(sig)
    sig = filtfilt(B,A,sig')';
else
    sig = complex(filtfilt(B,A,real(sig)')',filtfilt(B,A,imag(sig)')');
end
sig(nanInds) = nan;