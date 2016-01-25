function [posa posd fast w vel b] = fixPos(file,thresh)
%% take position info from whl file and preprocess it
% file is either a .whl file, or a tx4 array of positions
% posa = position array with '-1' removed
% posd = 1-d position on linear track, separating 2 directions of motion
% fast = boolean indicating if rat was faster than some threshold (useful
%   for removing periods when there is no motion and theta shrinks)
% w = trial number for linear track runs

if ~exist('thresh','var')
    thresh = .05;
end
bounds = [.05 .95];%[.15 .85];
if isstr(file)%strcmp(file,'char')
    %load('/media/work/hippocampus/KenjiData.mat');
    %whichDay = strcmp(file,Beh(:,4));
    %file = ['/media/Kenji_data/' Beh{whichDay,3} '/' Beh{whichDay,1} '/' file '/' file '.whl'];
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
posa = pos;
if nargout > 1
    nanInds = double(isnan(pos(:,1)) | isnan(pos(:,3)));
sz = round(size(pos,1)/2);nanInds(sz:end) = nanInds(sz:end)*2;
for i = 1:size(pos,2)
    pos(nanInds == 1,i) = pos(find(nanInds == 1, 1, 'last' )+1,i);
    pos(nanInds == 2,i) = pos(find(nanInds == 2,1)-1,i);
end
pos = bsxfun(@minus,pos,mean(pos));
vel = angVel(pos);
%vel1 = sqrt(sum(diff(pos(:,1:2)).^2,2));
%vel1 = filtLow(vel1,1250/32,1);
vel = filtLow(vel(:,1),1250/32,2);
[a,b,c] = svd(pos(:,1:2),'econ');pos = a;
for i = 1:2    
    pos(:,i) = pos(:,i) - prctile(pos(:,i),1);%min(pos(:,i));
    pos(:,i) = pos(:,i)/prctile(pos(:,i),99);%(max(pos(:,i)));
    pos(:,i) = max(0,min(pos(:,i),.9999));
end
posd = pos(:,1);
% %%FOR 1D TRACK
b = nan*ones(size(pos,1),1);
b(pos(:,1) < bounds(1)) = -1;b(pos(:,1) > bounds(2)) = 1;
nanInds = find(~isnan(b));
b = interp1(nanInds,b(nanInds),1:size(pos,1));
b = [0 diff(b)];
w = watershed(b==0);
w(w == 0) = w(find(w == 0)-1);
w = w-1; 
m = accumarray(w'+1,[0; diff(pos(:,1))],[],@mean);
posd(ismember(w+1,find(m>0))) = 1+posd(ismember(w+1,find(m>0)));
fast = vel > thresh*prctile(vel,99);%max(vel);
w = floor(w/2);%w = double(w-min(w) + 1); G removed this 1/10/2015
w(w == max(w)) = 0;%%G added this 1/10/2015
fast(w == 0) = 0;
end