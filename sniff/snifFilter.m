function a1 = snifFilter(a,f1,Fs,order)%,f2

%theta is 5-11Hz
%Fs = 1250;

if ~exist('order','var')
    order = 2;
end
%f1 = fdesign.bandpass('n,f3dB1,f3dB2',order,f1,f2,Fs);
%d1 = design(f1,'butter');
%a1 = filter(d1,flipud(filter(d1,flipud(a'))))';
%a1 = filter(d1,a')';
[d,n] = butter(order,f1/Fs*2,'high');
a1 = filtfilt(d,n,a')';