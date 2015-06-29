function [icasig,A,W,icaReshape,whichInds] = eyeCA(fn,trials,fullSess,pcs)

pad = 60; %1 sec
vis = 30;
%trials = 10:12;

load(fn,'all_cam2_ons','all_cam2_offs','lengths2','ims','sess'); %,'ims'

%ims = ims(:,:)';
%ims = bsxfun(@minus,ims,mean(ims,2));
if ~exist('pcs','var')
    pcs = 100;%size(v,1);%v,2);
end
%ims = ims + randn(size(ims))/30;
%figure;imagesc(squeeze(ims(500,:,:)));drawnow;
v = ims(:,:)';clear ims;
%v = v';

whichInds = logical(zeros(1,size(v,2)));
sessInds = zeros(1,size(v,2));
ctr = 1;
for i = 1:numel(trials)
    f = find(sess == trials(i));
    if ~fullSess
        for j = 1:numel(f)
            whichInds(lengths2(i)+(all_cam2_ons(f(j))-pad:all_cam2_offs(f(j))+pad)) = 1;
        end
    else
        whichInds(lengths2(i)+(all_cam2_ons(f(1))-pad:all_cam2_offs(f(end))+pad)) = 1;
    end
    for j = 1:numel(f)
        sessInds(lengths2(i)+((all_cam2_ons(f(j))-pad):(all_cam2_offs(f(j))+pad))) = ctr;
        ctr = ctr + 1;
    end
end
sessInds = sessInds(whichInds);
%[W,sphere,icasig,~,s] = runica(v(:,whichInds),'pca',pcs,'stop',5e-6);%,'extended',1
%W = W*sphere;
%A = pinv(W);

[icasig,A,W] = fastica(v(:,whichInds),'approach','symm','lastEig',pcs);%,'numOfIC',100,'g','tanh'
% The below is for using topographic / subspace ICA..
%   %p.model = 'tica';
%   %p.xdim = 10;
%   %p.ydim = 10;
%   %p.maptype = 'standard';
%   %p.neighborhood = 'ones3by3';
%   p.model = 'isa';
%   p.groupsize = 2;
%   p.groups = 50;
%   p.seed = 1;
%   %p.write = 5;
%   p.algorithm = 'gradient';
%   p.stepsize = 0.1;
%   p.epsi = 0.005;
%   [v,whiteningMatrix,dewhiteningMatrix] = whiten(v(:,whichInds));
%   whiteningMatrix = whiteningMatrix(end-(1:pcs)+1,:);
%   dewhiteningMatrix = dewhiteningMatrix(:,end-(1:pcs)+1);
%   v = v(end-(1:pcs)+1,:);
%   [icasig,A,W] = estimate(v, whiteningMatrix, dewhiteningMatrix, p );

icaReshape = zeros(size(icasig,1),max(sessInds),2*pad);
for i = 1:max(sessInds)
    icaReshape(:,i,:) = icasig(:,find(sessInds == i,1)+(0:size(icaReshape,3)-1));
end