function y = filtLow(sig,Fs,f,order)

if nargin < 4
    order = 8;
end
sig(isnan(sig)) = 0;
f1 = fdesign.lowpass('n,f3dB',order,f,Fs);
d1 = design(f1,'butter');
%a1 = filter(d1,flipud(filter(d1,flipud(a'))))';
[B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
if isreal(sig)
    y = filtfilt(B,A,sig')';
else
    y = complex(filtfilt(B,A,real(sig)')',filtfilt(B,A,imag(sig)')');
end