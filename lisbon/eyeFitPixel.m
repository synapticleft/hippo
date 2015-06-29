function [cc,mse,kern] = eyeFitPixel(fn)
% Do a regression to fit each pixel in a movie, based on pupil fits
% (position and size)
ridges = 10.^(-8:-2);
load(fn,'pupils','trialInds','ims','lengths2');

temp = zeros(size(pupils,1),numel(lengths2)-1);
for i = 1:numel(lengths2)-1
    temp(lengths2(i)+1:lengths2(i+1),i) = 1;
end

nanInds = ~isnan(pupils(:,1));

ims = zscore(ims(nanInds,:));
pupils = [pupils(:,1:5) temp];

pupils = pupils(nanInds,:);
pupils = taylorExpand(pupils,1:5);%
pupils(:,end+1) = 1;

% w = ims'/pupils';
% imsHat = w*pupils';
% cc = diag(imsHat*ims);
% cc = cc/sum(nanInds);
[cc,mse,kern] = ridgeCross(ims,pupils,3,ridges);
figure;imagesc(reshape(cc,[30 40]));


function y = taylorExpand(y,inds)
if ~exist('inds','var')
    inds = 1:size(y,2);
end
for i = inds
    y = [y bsxfun(@times,y(:,inds),y(:,i))];
end