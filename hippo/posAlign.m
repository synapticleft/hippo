function allRuns = posAlign(pos,anchor)

margin = [-250 100];
[~,p1,~,tr] = fixPos(pos);
%p1(p1 > 1) = 2-p1(p1 > 1);
allRuns = zeros(max(tr),range(margin)+1);
figure;
for i = 1:max(tr)
    pt = p1(tr == i);
    [~,ind] = min(abs(pt-anchor));
    try
    %plot(pt(ind+(margin(1):margin(2))));hold all;
    %allRuns(i,:) = pt(ind+(margin(1):margin(2)));
    temp = p1(find(tr == i,1) + ind+(margin(1):margin(2)));
    temp(temp > 1) = 2-temp(temp>1);
    plot(temp);hold all;
    allRuns(i,:) = temp;
    catch
        i
    end
end