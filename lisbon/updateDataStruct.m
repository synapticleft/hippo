function data1 = updateDataStruct

matFiles = {'cam2_7-9Dec16.mat','cam2_10-12Dec16.mat','cam2_13-15Dec16.mat'};
sess = {4:6,7:9,10:12};
datFile = 'data_s1_corr2.mat';

load(datFile,'data','all_cam2_ons','all_cam2_offs');
data1{1,:} = {'PupilX Fit','PupilY Fit','PupilMajor Fit','PupilMinorFit','PupilTheta Fit',...
    'PupilX Fit HP','PupilY Fit HP','PupilMajor Fit HP','PupilMinorFit HP','PupilTheta Fit HP'};

Fs = 30;
highPass = .01;
ridge = .01;
for i = 1:numel(matFiles)
    load(matFiles{i},'pupils','ims');
    ims = bsxfun(@minus,ims,mean(ims));
    nanInds = ~isnan(pupils(:,1));
    xCov = ims(nanInds,:)'*ims(nanInds,:);
    xCov = xCov + ridge*diag(diag(xCov));
    w = pupils(nanInds,1:5)'*ims(nanInds,:)/xCov';
    pupilsHat = w*ims(:,:)';
    pupilsHath = filtHigh(pupilsHat,Fs,highPass);
    for j = 1:numel(sess{i})
        f = find([data{2:end,2}] == sess{i}(j));
        for k = 1:numel(f)
            inds = all_cam2_ons(f(k))-60:all_cam2_offs(f(k))+58;
            for l = 1:size(pupilsHat,1)
                data1{f(k)+1,l} = pupilsHat(l,inds);
                data1{f(k)+1,l+size(pupilsHat,1)} = pupilsHath(l,inds);
            end
            data1{f(k)+1,2*size(pupilsHat,1)+1} = interp1(data{f(k)+1,22},data{f(k)+1,15},data{f(k)+1,28},'pchip');
            data1{f(k)+1,2*size(pupilsHat,1)+2} = nan*ones(size(data{f(k)+1,28}));
            fnan = find(~isnan(data{f(k)+1,6}));
            temp = interp1(data{f(k)+1,22}(fnan),data{f(k)+1,6}(fnan),data{f(k)+1,28}(find(data{f(k)+1,28} == 0):...
                find(data{f(k)+1,28} <= data{f(k)+1,22}(fnan(end)),1,'last')));
            data1{f(k)+1,2*size(pupilsHat,1)+2}(find(data{f(k)+1,28} == 0)+(0:numel(temp)-1)) = temp;
        end
    end
end