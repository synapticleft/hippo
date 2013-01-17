function y = filtNotch(sig,Fs,f,order,q)

if nargin < 3
    f = 60;
    order = 4;
    q = 100;
end
nanInds = isnan(sig);
sig(nanInds) = 0;
f1 = fdesign.notch('n,f0,q',order,f,q,Fs);
d1 = design(f1,'butter');
%a1 = filter(d1,flipud(filter(d1,flipud(a'))))';
[B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
if isreal(sig)
    y = filtfilt(B,A,sig')';
else
    y = complex(filtfilt(B,A,real(sig)')',filtfilt(B,A,imag(sig)')');
end
y(nanInds) = nan;