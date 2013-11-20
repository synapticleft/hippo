function Xf = morFilter(X,Fc,Fs,norm,win)

%X = n x t matrix, where n = # channels, t = # timesteps
%Fc = center frequency
%Fs = sampling rate
%norm = normalize output of each channel to 1 (default 0)
%win = window for filter (default 1 s)
%Xf = complex-valued filtered output

if ~exist('norm','var')
    norm = 0;
end
if ~exist('win','var')
    win = 1;
end
Fb = 1/500; %This parameter controls bandwidth, may need to be changed for your purposes
%Fc = 8.4 for hippocampus
%Fb = 1/500, win = 1 for hippocampal LFP
%Fb = 1/5000, Fc = 160, win = 5 for ripples
[psi,x] = cmorwavf(-win,win,Fs*2*win,Fb,Fc);
%figure;plot(x,real(psi));hold all;plot(x,imag(psi));plot(x,abs(psi));return
Xf = zeros(size(X));
for i = 1:size(X,1)
    Xf(i,:) = fliplr(filter(conj(psi),1,fliplr(filter(psi,1,X(i,:)))));
    if (norm)
        Xf(i,:) = Xf(i,:)/std(Xf(i,:));
    end
end