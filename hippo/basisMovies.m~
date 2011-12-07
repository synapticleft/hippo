function basisMovies(file,gridDim,ratio)
%makes movies of 2-d basis functions running in time

curDir = ['/media/work/hippocampus/'];
phi = hdf5read([curDir file],'phi');%h5varget
order = hdf5read([curDir file],'order');
variance = hdf5read([curDir file],'variance');
[junk sorted] = sort(variance,'descend');
%sorted = order + 1;

while 1
    for i = 1:size(phi,1)
        tic;
        showGrid(squeeze(phi(i,sorted,:))',gridDim,ratio);
        pause(.1-toc);
    end
end

