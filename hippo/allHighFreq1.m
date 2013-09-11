function [A,tes,td] = allHighFreq1(file,pos)

X = getData(file,1,1:256);
[size(X) size(pos)]
X = filtHigh(X,1250,100,4);
[A,tes,td] = runHighFreq(X,pos);
save('512HighFreq.mat','A','tes','td');