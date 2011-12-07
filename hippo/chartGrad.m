function [temp absCourse] = chartGrad(data,bins)

if numel(bins)  == 1
    bins(2) = bins(1);
end
[u s v] = svds(data,2);
data1 = u*s*v';

[~,absRank] = sort(abs(v(:,1)));
%cc = corr(circ_dist(circ_mean(angle(u(:,1))),angle(u(:,1))),circ_dist(angle(data1),angle(u(:,1)*v(:,1)')),'type','kendal');%'spearman'
angDev = circ_dist(angle(data),angle(u(:,1)*v(:,1)'));
cc = myCircCorr(angle(u(:,1)),angDev);
[~,gradRank] = sort(cc);
% figure;
% lims{1} = linspace(min(angle(u(:,1))),max(angle(u(:,1))),50);%
% lims{2} = linspace(-1,1,20)/3;
% for i = 1:bins(1)
%     for j = 1:bins(2)
%         r1 = max(1,ceil((i-1)/bins(1)*numel(absRank))):ceil(i/bins(1)*numel(absRank));
%         r2 = max(1,ceil((j-1)/bins(2)*numel(absRank))):ceil(j/bins(2)*numel(absRank));
%         inds = intersect(absRank(r1),gradRank(r2));%absRank == i & gradRank == j;
%         subplot(bins(1),bins(2),(i-1)*bins(2)+j);
%         temp = angDev(:,inds);
%         imagesc(lims{1},lims{2},log(hist3([repmat(angle(u(:,1)),[numel(inds) 1]) temp(:)],lims))');% circ_dist(angle(data1(:,inds)),angle(u(:,1)*v(inds,1)'))']));drawnow;
%         allDat(i,j,:) = circ_mean(temp');
%     end
% end
% figure;plot(angle(u(:,1)),squeeze(mean(allDat,1))','.');
% figure;plot(angle(u(:,1)),squeeze(mean(allDat,2))','.');

angCourse = sum(bsxfun(@times,angDev,circ_dist(circ_mean(angle(u(:,1))),angle(u(:,1)))));
angAbsCourse = sum(bsxfun(@times,bsxfun(@minus,abs(data),mean(abs(data))),circ_dist(circ_mean(angle(u(:,1))),angle(u(:,1)))));
%stdCourse = circ_std(angle(data1));
figure;plot(angCourse/std(angCourse));hold all;plot(abs(v(:,1))/std(abs(v(:,1))));plot(angAbsCourse);
absCourse = abs(v(:,1));
temp = absCourse.*angCourse';
plot(temp/std(temp));%plot(stdCourse,'r');
figure;scatter(angCourse,absCourse);
figure;plot(xcov(temp,200,'coeff'));hold all;plot(xcov(abs(v(:,1)),200,'coeff'));
%figure;scatter(angCourse,log(stdCourse));