function getConFit(dat,thresh)

magx = 10;
for i = 1:size(dat,1)
    for j= 1:size(dat,2)
        coords(size(dat,1)*(i-1) + j,:) = [i,magx*j,dat(i,j)];
    end
end
mask = coords(:,3) > thresh;
%figure;scatter3(coords(:,1),coords(:,2),coords(:,3));
figure;hold on;scatter3(coords(mask,1),coords(mask,2),coords(mask,3),'r');
[a b phi r e f] = lscone(coords,[-30 30 30]',[0 0 1]',pi/4,50,1,1,mask);
%figure;hold on;
a
thetas = 0:.3:(2*pi);
for i = 0:10
    for j = 1:length(thetas)
        scatter3(a(1) + (tan(phi)*cos(thetas(j))+b(1))*i,a(2) + (tan(phi)*sin(thetas(j))+b(2))*i,(a(3)-b(3)*i)+4,'b');
    end
end
    

%

        
        