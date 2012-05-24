%% get some information
data_dir = '/media/work/lgn/';
allFiles = dir(data_dir);
data_file = allFiles(fileNum+2).name;
%data_file = '040218E_combined.3.h5';
%data_file = '050204M_mbmonkcc000.3.h5';
%data_file = '050921J_jecatcc000.3.h5';
%data_file = '050921L_lmcatcc000.3.h5';
%data_file = '040218D_combined.3.h5';
fname = [data_dir data_file];
hinfo = hdf5info(fname);

%% list groups
%hinfo.GroupHierarchy.Groups.Name

%% list datasets in first group
%hinfo.GroupHierarchy.Groups(1).Datasets.Name

%% read raw data

signals = hdf5read(fname, '/data/signals');
%size(signals) % (29473, 20)

%% read spike and epsp train

spktrain = hdf5read(fname, '/data/spktrain');
psptrain = hdf5read(fname, '/data/psptrain');
signals = hdf5read(fname, '/data/signals');

%% read frequency and band width

fm = hdf5read(fname, '/data/f');
sf = hdf5read(fname, '/data/sf');

%% read analytic signal
% this stalls matlab fore some reason
%asigs = hdf5read(fname, '/data/asigs');
asigs = h5varget(fname, '/data/asigs');
asigs = complex(asigs.r,asigs.i);