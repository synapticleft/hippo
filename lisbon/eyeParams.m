function [ims,shifts] = eyeParams(dn) 
% This function is meant to find distortions in the eye, and track them
% over time for the purposes of feature extraction

%considerations
% - align to previous frame, or to first in sequence? previous frame seems
% more intuitive?
% - align only frames within trial, or all frames in the session?
% - downsample? 4 seems reasonable..?
% - demons - # pyramids, iterations, smoothness?

%compare the numbers we get for downsampling, aligning to prev vs first
prev = 1;
scale = 4;
sz = [480 640]/scale;
cd(['/media/gaba/Elements/work/lisbon/Roberto Data/' dn]);
file = dir('d*.mat');
load(file(1).name,'data');
file = dir('E*.mat');
load(file(1).name);
file = dir('*.avi');
ims = [];
sess = [data{2:end,2}];
ctr = 1;
lens = all_cam2_offs-all_cam2_ons;
ims = zeros(sum(lens(sess == numel(file))),sz(1),sz(2));
for i = numel(file)%1%:
    currSess = find(sess == i);
    m = VideoReader(file(i).name);
    for j = currSess%(1)    
       for k = all_cam2_ons(j):all_cam2_offs(j)-1
           m.CurrentTime = (k+.5)/m.FrameRate;
           %ims(ctr,:,:) = imresize(rgb2gray(readFrame(m)),1/scale);
           shifts(ctr) = j;
           ctr = ctr + 1;
       end
%         %%ALIGNMENT CODE
%         m.CurrentTime = (all_cam2_ons(j)-.5)/m.FrameRate;
%         oldFrame = imresize(rgb2gray(readFrame(m)),1/scale);
%         ims(1,:,:) = oldFrame;
%         for k = all_cam2_ons(j):all_cam2_offs(j)-1
%             tic;
%             m.CurrentTime = (k+.5)/m.FrameRate;
%             newFrame = imresize(rgb2gray(readFrame(m)),1/scale);
%             newFrame = imhistmatch(newFrame,oldFrame);
%             currInd = k-all_cam2_ons(j)+2;
%             [shifts(currInd,:,:,:)] = imregdemons(newFrame,oldFrame);%,ims(currInd,:,:)
%             ims(currInd,:,:) = newFrame;
%             if prev
%                 oldFrame = newFrame;
%             end
%             %subplot(211);imshowpair(oldFrame,newFrame);axis image;
%             %subplot(212);imshowpair(oldFrame,squeeze(ims(currInd,:,:)));axis image;drawnow;
%             %toc
%        end
%         %% END ALIGNMENT CODE
    end
end