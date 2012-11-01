function [Xf Xf1] = demodFilter(X,Fc,Fs,norm)
if ~exist('norm','var')
    norm = 0;
end
win = 1;%/5;
Fb = 1/500;
%Fc = 8.4 for hippocampus
%Fb = 1/500, win = 1 for hippocampal LFP
%Fb = 1/5000, Fc = 160, win = 5 for ripples
[psi,~] = cmorwavf(-win,win,Fs*2*win,Fb,Fc);
[psi1,~] = cmorwavf(-win,win,Fs*2*win,Fb*20,Fc);
%figure;subplot(211);plot(real(psi));hold all;plot(imag(psi),'r');
%subplot(212);plot(real(psi1));hold all;plot(imag(psi1),'r');
Xf = zeros(size(X));
for i = 1:size(X,1)
%    Xf(i,:) = filter(psi,1,X(i,:));
    Xf(i,:) = fliplr(filter(conj(psi),1,fliplr(filter(psi,1,X(i,:)))));
%    Xf1(i,:) = fliplr(filter(conj(psi1),1,fliplr(filter(psi1,1,X(i,:)))));
%    if (norm)
%        Xf(i,:) = Xf(i,:)/std(Xf(i,:));
%    end
end
[u s v] = svds(Xf,1);
v1 = fliplr(filter(conj(psi1),1,fliplr(filter(psi1,1,v'))));
figure;plot(circ_diff(v(1:1000)));hold on;plot(circ_diff(v1(1:1000)),'r');
Xf1 = Xf.*exp(1i*-angle(u*v1));
%Xf1 = bsxfun(@times,Xf,exp(1i*-angle(u*v1)));