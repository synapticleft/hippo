function pupils = eyeExtract(fns,outputFile)
%extract details from eye file triplets and save into a .mat file

matData = 'data_s1_corr2.mat';

aviFiles = dir('*.avi');
csvFiles = dir('*.csv');

load(matData,'data','all_cam2_ons','all_cam2_offs');

sess = [data{2:end,2}];
lengths2 = 0;
ims = [];
pupils = [];
for i = 1:numel(fns)
   m = VideoReader(aviFiles(fns(i)).name);
    c = csvread(csvFiles(fns(i)).name);
    lengths2 = [lengths2 size(c,1)];
    readFrame(m);
    imst = zeros(size(c,1),30,40);
    for j = 1:size(c,1)
        imst(j,:,:) = imresize(rgb2gray(double(readFrame(m))/255),1/16);
    end
    ims = [ims;imst];
    pupils = [pupils; c];
end

pupils(pupils == 0) = nan;
pupils = bsxfun(@minus,pupils,nanmean(pupils));
pupils = bsxfun(@rdivide,pupils,nanstd(pupils));
lengths2 = cumsum(lengths2);
save(outputFile,'sess','ims','lengths2','pupils');