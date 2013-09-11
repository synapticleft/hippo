function a1 = hipFilter(a,f1,f2,Fs,order)
%% old way of filtering hippo data

%theta is 5-11Hz
%Fs = 1250;

if nargin < 5
    order = 8;
end
f1 = fdesign.bandpass('n,f3dB1,f3dB2',order,f1,f2,Fs);
d1 = design(f1,'butter');
%a1 = filter(d1,flipud(filter(d1,flipud(a'))))';
a1 = filter(d1,a')';
%a1 = filtfilt(d1.sosMatrix,d1.ScaleValues,a')';

% a1 = zeros(size(a));
% 
% for i = 1:size(a,1)
%     a1(i,:) = filter(d1,squeeze(a(i,:)));
% end
