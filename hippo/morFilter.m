function Xf = morFilter(X,Fc,Fs)%,win

%X = n x t matrix, where n = # channels, t = # timesteps
%Fc = center frequency
%Fs = sampling rate
%norm = normalize output of each channel to 1 (default 0)
%win = window for filter (default 1 s)
%Xf = complex-valued filtered output

%if ~exist('norm','var')
%    norm = 0;
%end
%if ~exist('win','var')
%    win = 2;
%end
%Fb = 1/(2*pi*Fc.^2);%500; %This parameter controls bandwidth, may need to be changed for your purposes
%Fc = 8.4 for hippocampus
%Fb = 1/500, win = 1 for hippocampal LFP
%Fb = 1/5000, Fc = 160, win = 5 for ripples
%psi = mf_cmorlet(Fc,Fs);
Fb = 1/500;
win = 1;
[psi,x] = cmorwavf(-win,win,Fs*2*win,Fb,Fc);
inds = abs(psi) > .001*max(abs(psi));
psi = psi(inds);x = x(inds);
psi = psi/sqrt(sqrt(sum(abs(psi).^2)));
%Xf = fconv(X,psi);
 Xf = zeros(size(X));
 for i = 1:size(X,1)
     Xf(i,:) = fliplr(filter(conj(psi),1,fliplr(filter(psi,1,X(i,:)))));
 end

function y = mf_cmorlet(f0,fs)
ts = 1/fs; %sample period, second
SD_f = f0/5;
SD_t = 1/(2*pi*SD_f);
t = [0:ts:4*SD_t];
t = [-t(end:-1:2) t];
y = exp( -t.^2/(2*SD_t^2) ).* exp( 1i*2*pi*f0*t );
y = y/sqrt(SD_t*sqrt(pi));%(sqrt(0.5*sum( real(y).^2 + imag(y).^2 )));%