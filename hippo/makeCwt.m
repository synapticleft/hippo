function makeCwt(file,elecs,pos)
%% make a CWT for either spikes or LFP

lfp = 0;
dec = 32;
%elecs = 33:34;%64;
levels = 1:.5:8;
someShanks = 1:8;bin = 1/1.25;

if nargin < 3
    pos = [file '.whl'];
end
[~,~,~,runNum] = fixPos(pos);

runNum = round(interp(runNum,dec));

Xc = single(zeros(numel(levels),numel(runNum)));
if lfp
    fid = fopen([file 'Cwt.dat'],'a');
else
    fid = fopen([file 'CwtSpk.dat'],'a');
    [spikeTime,spikeID,~,info] = LoadCluRes(file,someShanks);
    spikeTime = ceil(spikeTime/info.SampleRate*1000/bin);
    [~, cellInfo] = hipSpikes(file,bin,'',1);
    spikeTime(ismember(spikeID,find(cellInfo.type < 1))) = [];
    spikeID(ismember(spikeID,find(cellInfo.type < 1))) = [];
    %allIDs = unique(spikeID);
    h = hist(spikeID,1:max(spikeID));
    elecs = find(h > 500);
end

for i = elecs
    if lfp 
    X = getData(file,1,i);
    if numel(runNum) < numel(X)
        X = X(1:numel(runNum));
    else
        runNum = runNum(1:numel(X));
    end
    else
        X = zeros(size(runNum));
        X(spikeTime(spikeID == i)) = 1;
    end
    for j = 1:max(runNum)
        Xc(:,runNum == j) = single(cwt(X(runNum == j),2.^levels,'cmor1-1'));
    end
    fwrite(fid,[real(Xc);imag(Xc)],'single');
    i
end
fclose(fid);
