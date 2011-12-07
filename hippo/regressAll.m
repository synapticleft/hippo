function regressAll(pcaDat,lfpDat,posDat,inds)

if size(inds) == 1
    inds = 10000+(1:inds);
end
if size(lfpDat,1) > size(lfpDat,2)
    lfpDat = lfpDat.';
end
if size(pcaDat,1) > size(pcaDat,2)
    pcaDat = pcaDat.';
end
if size(posDat,1) > size(posDat,2)
    posDat = posDat.';
end
lfpDat = lfpDat(:,inds);pcaDat = pcaDat(:,inds);posDat = posDat(:,inds);
relPhase = diff(angle(pcaDat(1:2,:)));
instPhase = [diff(angle(pcaDat),1,2) zeros(size(pcaDat,1),1)];
input = [pcaDat;relPhase;exp(1i*relPhase);abs(pcaDat(1,:)).*exp(1i*relPhase);...
    abs(pcaDat(2,:)).*exp(1i*relPhase);circ_std(angle(lfpDat));instPhase];
%input = (pcaDat(1,:));
input = posDat;
output = abs(pcaDat(1:2,:));
%output = (pcaDat(2,:));
%output = [pimpOut(posDat,1); pimpOut(posDat,2)];
%output = output(9,:);
%output = filtLow(output,1250/32,6);
%input = [input;input.*conj(input)];
%input = output(1:10,:);
%output = output(13,:);
offset = 0;
[c m] = splitToep1(output,input,50,1250/32,3,size(input,2)-offset,offset,1000);
if numel(c) == 1
    c
else
    figure;plot(abs(c))
end

function p = pimpOut(data,nDiff)
data(data == -1) = nan;
orient = angle(complex(diff(data([1 3],:)),diff(data([2 4],:))));
if nDiff 
    data = [diff(data,nDiff,2) zeros(size(data,1),nDiff)];
end
data(isnan(data)) = 0;orient(isnan(orient)) =0;
cData = complex(data([1 3],:),data([2 4],:));
relC = bsxfun(@times,cData,exp(-1i*orient));
%figure;hist(angle(relC(1,:)),100)
%circ_std(angle(relC'))
p = [data;cData;exp(1i*angle(cData));abs(cData);relC;real(relC);imag(relC);diff(cData);abs(diff(cData))];