function makeCwt(file,elecs,pos)

dec = 32;
%elecs = 33:34;%64;
levels = 1:.5:8;

if nargin < 3
    pos = [file '.whl'];
end
[~,~,~,runNum] = fixPos(pos);

runNum = round(interp(runNum,dec));

Xc = zeros(numel(levels),numel(runNum));
fid = fopen([file 'Cwt.dat'],'a');
for i = elecs
    X = getData(file,1,i);
    if numel(runNum) < numel(X)
        X = X(1:numel(runNum));
    else
        runNum = runNum(1:numel(X));
    end
    for j = 1:max(runNum)
        Xc(:,runNum == j) = cwt(X(runNum == j),2.^levels,'cmor1-1');
    end
    fwrite(fid,[real(Xc);imag(Xc)],'single');
    i
end
fclose(fid);