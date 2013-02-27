function [dx dy] = angGradient(in)
dx = ones(size(in));dy = ones(size(in));
dx(:,2:end) = in(:,2:end)./in(:,1:end-1);
dx(:,1:end-1) = dx(:,1:end-1).*dx(:,2:end);
dx = angle(dx);
dy(2:end,:) = in(2:end,:)./in(1:end-1,:);
dy(1:end-1,:) = dy(1:end-1,:).*dy(2:end,:);
dy = angle(dy);