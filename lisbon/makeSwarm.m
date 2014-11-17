function makeSwarm(fn,inds,trs,whichDiffs)

if ~exist('whichDiffs','var')
    whichDiffs = [-4 -2 -.5 .5 2 4];
end
if exist('trs','var')
    [allData a] = preProcessRoberto(fn,inds,0,0,trs);
else
    [allData a] = preProcessRoberto(fn,inds,0,0);
end

figure;
for i = 40:180
for j = 1:numel(whichDiffs)
subplot(2,numel(whichDiffs)/2,j);
f = find(a(:,6) == whichDiffs(j) & a(:,3) == -1);
scatter(allData(f,i,1),allData(f,i,2),'b','filled');hold on;
scatter(allData(f,i,3),allData(f,i,4),'r','filled');
%plot(squeeze(allData(f,i,[1 3]))',squeeze(allData(f,i,[2 4]))','b');
hold off;
xlim([-3 3]);ylim([-2 2]);
end
drawnow;
m(i-39) = getframe(gcf);
end
movie2avi(m,'swarm1.avi');