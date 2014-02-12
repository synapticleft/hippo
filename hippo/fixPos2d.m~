function [pos] = fixPos2d(file)
%% take position info from whl file and preprocess it
% file is either a .whl file, or a tx4 array of positions
% posa = position array with '-1' removed
% posd = 1-d position on linear track, separating 2 directions of motion
% fast = boolean indicating if rat was faster than some threshold (useful
%   for removing periods when there is no motion and theta shrinks)
% w = trial number for linear track runs

thresh = .05;
if isstr(file)
    pos = importdata(file);
else
    pos = file;
end

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