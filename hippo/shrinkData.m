function shrinkData(fileIn,suffix)
%% take dat file, decimate variables, and save them as variables in a single mat file.

fileOut = fileIn;
decs = [ones(1,3) 2*ones(1,5) 8*ones(1,7)];
dec = 32;

levels = 1:.5:8;
d = ['/media/work/hippocampus/' fileIn '/'];
m = memmapfile([[d fileIn] suffix '.dat']);
dims = [30 4947136];%1560096];
dims(3) = numel(m.Data)/prod(dims)/4;

m = memmapfile([[d fileIn] suffix '.dat'],'Format',{'single' dims 'X'});
[~,pos,isFast,~] = fixPos([fileIn '.whl']);
for i = numel(levels):-1:1
    %X1 = permute(m.Data(1).X([i i+numel(levels)],:,:),[2 1 3]);
    isFast1 = logical(round(interp([double(isFast); 0],dec/decs(i))));
    X = single(zeros(dims(3),sum(isFast1)));
    for j = 1:dims(3)
        X1 = squeeze(m.Data(1).X([i i+numel(levels)],:,j))';
        temp = decimate(squeeze(double(complex(X1(:,1),X1(:,2)))),decs(i));
        X(j,:) = single(temp(isFast1));
    end
    save([fileOut suffix num2str(i) '.mat'],'X');
    i
end