function [v h] = angVel(pos)
%% find velocity of rat with respect to head orientation

pos(pos == -1) = nan;
or = pos(:,1:2) - pos(:,3:4);
or(end,:) = [];
pos = diff(pos);
h = hypot(or(:,1),or(:,2));
or = bsxfun(@rdivide,or,h);

v(:,1) = sum(pos(:,1:2).*or,2);
v(:,3) = sum(pos(:,3:4).*or,2);
or(:,2) = -or(:,2);or = or(:,[2 1]);
v(:,2) = sum(pos(:,1:2).*or,2);
v(:,4) = sum(pos(:,3:4).*or,2);

for i = 1:size(v,2)
    nanInds = find(~isnan(v(:,i)));
    v(:,i) = interp1(nanInds,v(nanInds,i),1:size(v,1));
end
%or(:,2) = -or(:,2);or = or(:,[2 1]);

%figure;plot(v(:,[1 2]));