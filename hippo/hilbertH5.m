function hilbertH5(input,output,numChans)

h5filecreate(output);
[~,~,nSamples] = LoadBinary(input,1,[],[],[],[],[1 2]);
%nSamples = 10000;
h5datacreate(output,'/hReal','type','int16','size',[numChans nSamples]);
h5datacreate(output,'/hImag','type','single','size',[numChans nSamples]);
for i = 1:numChans
a = LoadBinary(input,i,[],[],[],'int16',[1 nSamples]);
h = hilbert(single(a));
h5varput(output,'/hReal',[i-1 0],[1 nSamples],real(h));
h5varput(output,'/hImag',[i-1 0],[1 nSamples],imag(h));
i
end