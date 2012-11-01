function v2 = timeScale(v,v1)

av = angle(v(:,1));
v1(:,1) = v1(:,1).*exp(1i*-1);
av1 = angle(v1(:,1))/2;
%av2 = angle(exp(1i*(angle(v1(:,1))/2+pi)));
avs = [angle(v1(:,1))/2 angle(exp(1i*(angle(v1(:,1))/2+pi)))];
%smaller = abs(circ_dist(circ_dist(av,avs(:,1)),-.5)) < abs(circ_dist(circ_dist(av,av2),-.5));
smaller = abs(circ_dist(av,avs(:,1))) < abs(circ_dist(av,avs(:,2)));
avn(smaller) = avs(smaller,1);avn(~smaller) = avs(~smaller,2);
v2 = abs(v1(:,1)).*exp(1i*avn');%min(abs(circ_dist(circ_dist(av,av1),-.5)),abs(circ_dist(circ_dist(av,av2),-.5))));
hist(circ_diff(v2),100);
%figure;hist(circ_dist(av,angle(v2)),100);
%sPlot(([v(:,1) v1(:,1) v2]'));
%figure;imagesc(hist3([circ_dist(av,av1) abs(v1(:,1))],[100 100]));