function Xf = morFilter(X,Fc,Fs,norm)
if ~exist('norm','var')
    norm = 1;
end
win = 1;%/5;
Fb = 1/500;
%Fc = 8.4 for hippocampus
%Fb = 1/500, win = 1 for hippocampal LFP
%Fb = 1/5000, Fc = 160, win = 5 for ripples
[psi,x] = cmorwavf(-win,win,Fs*2*win,Fb,Fc);
Xf = zeros(size(X));
for i = 1:size(X,1)
    Xf(i,:) = fliplr(filter(conj(psi),1,fliplr(filter(psi,1,X(i,:)))));
    if (norm)
        Xf(i,:) = Xf(i,:)/std(Xf(i,:));
    end
end