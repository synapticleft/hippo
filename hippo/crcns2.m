function crcns2(traceT,dims,ratio,im,fileName)
% this function takes in dynamical data (traceT), dimensions of array
% (dims), relative spacing in horizontal and vertical directions (ratio),
% and a background image (im), to produce a movie. This movie can be saved
% to a file.

ds = 5;
subplot('position',[0 0 1 1]);
figure('Name','Traveling LFP waves in rat hippocampus');
if numel(dims) > 2
    probes = dims-min(dims(:))+1;%trace = trace(dims-min(dims(:)) + 1,:);
    dims = size(dims);
end
%numPast = 500;
%trange = std(traceT(:))*[-2 2];
temp = traceT(:,1);temp = temp(probes);
s = surf((1:dims(1))*ratio(1),(1:dims(2))*ratio(2),rot90(rot90(rot90(temp))),'edgecolor','none');
hold on;
surf([1 dims(1)]*ratio(1),[1 dims(2)]*ratio(2),repmat(-14, [2 2]),im,'facecolor','texture');
set(gca,'xlim',[1 dims(1)]*ratio(1),'ylim',[1 dims(2)]*ratio(2),'zlim',[-14 14]);
ylabel ('posterior-anterior');xlabel('ventral-dorsal');zlabel('voltage');
set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[],'fontsize',16);
axesLabelsAlign3D;
for i = 1:size(traceT,2)%startInd + 500%
    temp = traceT(:,i);temp = temp(probes);
    temp = imfilter(temp,fspecial('gaussian',5,1));
        set(s,'ZData',rot90(rot90(rot90(temp))));
        title([num2str((i-1)*ds/1250,2) ' s']);
        drawnow;
        %pause(.03);
        if exist('fileName','var')
            m(i) = getframe(gcf);
        end
end
if exist('fileName','var')
    movie2avi(m,fileName);
end

