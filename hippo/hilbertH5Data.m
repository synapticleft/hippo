function hilbertH5Data(Xf,output,numChans)

[numChans nSamples] = size(Xf);
%Xf = morFilter(X,8.4,1250/32);
h5filecreate(output);
h5datacreate(output,'/hReal','type','double','size',[numChans nSamples]);
h5datacreate(output,'/hImag','type','double','size',[numChans nSamples]);
for i = 1:numChans
    h5varput(output,'/hReal',[i-1 0],[1 nSamples],real(Xf(i,:)));
    h5varput(output,'/hImag',[i-1 0],[1 nSamples],imag(Xf(i,:)));
end