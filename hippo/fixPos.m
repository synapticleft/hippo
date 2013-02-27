function pos = fixPos(file)

pos = importdata(file);

pos(pos == -1) = nan;
reject = 0;
for i = 1:size(pos,2)
    reject = reject | min([0; diff(pos(:,i))],flipud([0; diff(flipud(pos(:,i)))])) < -20;
end
pos(reject,:) = nan;
for i = 1:size(pos,2)
    nanInds = find(~isnan(pos(:,i)));
    pos(:,i) = interp1(nanInds,pos(nanInds,i),1:size(pos,1));
end