function  realH5(input,output,numChans)

h5filecreate(output);
[~,~,nSamples] = LoadBinary(input,1,[],[],[],[],[1 2]);
%nSamples = 10000;
h5datacreate(output,'/hReal','type','int16','size',[numChans nSamples]);
for i = 1:numChans
a = LoadBinary(input,i,[],[],[],'int16',[1 nSamples]);
h5varput(output,'/hReal',[i-1 0],[1 nSamples],a);
i
end