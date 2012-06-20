function X = load_datachunk(m,p)
if strcmp(p.data_type,'lfp')
%    decFac = 1;
    info = hdf5info([p.data_root p.data_file]);
    nSamples = info.GroupHierarchy.Datasets(1).Dims;
     padding = 0;%20000;
    rind = padding + ceil(rand*(nSamples(2) - padding*2 - p.imszt));%*decFac
%     if decFac > 1
%          chunk1 = complex(double(h5varget([p.data_root 'hippo.h5'],'/hReal',[0 rind-1],[nSamples(1) decFac*p.imszt])),...
%              double(h5varget([p.data_root 'hippo.h5'],'/hImag',[0 rind-1],[nSamples(1) decFac*p.imszt])));
%          for i = 1:size(chunk1,1)
%              X(i,:) = decimate(chunk1(i,:),decFac);
%          end
%     else
%         X = double(h5varget([p.data_root p.data_file],p.data_field,[0 rind-1],[nSamples(1) p.imszt]));%'hippo.h5'
         X = complex(double(h5varget([p.data_root p.data_file],'/hReal',[0 rind-1],[nSamples(1) p.imszt])),...
            double(h5varget([p.data_root p.data_file],'/hImag',[0 rind-1],[nSamples(1) p.imszt])));
%     end
elseif strcmp(p.data_type,'sim')
    X = real(makeWavesMor(m,[m.patch_sz p.imszt]));
end