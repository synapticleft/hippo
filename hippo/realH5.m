function  realH5(input,output,numChans)

h5filecreate(output);
if exist('numChans','var')
    [~,~,nSamples,numChans] = LoadBinary(input,1,numChans,[],[],[],[1 2]);
else
    [~,~,nSamples,numChans] = LoadBinary(input,1,[],[],[],[],[1 2]);
end
[nSamples numChans]
%numChans = 96;
step = 400000;
h5datacreate(output,'/hReal','type','int16','size',[numChans nSamples]);
for i = 1:step:nSamples
    numElems = min(step,nSamples-i+1);
    a = LoadBinary(input,1:numChans,numChans,[],[],'int16',i+[1 numElems]-1);
    h5varput(output,'/hReal',[0 i-1],[numChans numElems],a);
    i
end
% 
% h5filecreate(output);
% [~,~,nSamples] = LoadBinary(input,1,[],[],[],[],[1 2]);
% %nSamples = 10000;
% h5datacreate(output,'/hReal','type','int16','size',[numChans nSamples]);
% for i = 1:numChans
% a = LoadBinary(input,i,[],[],[],'int16',[1 nSamples]);
% h5varput(output,'/hReal',[i-1 0],[1 nSamples],a);
% i
% end