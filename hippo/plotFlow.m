function plotFlow(u, v, imgOriginal)
% u and v are the horizontal and vertical optical flow vectors,
% respectively. imgOriginal, if supplied, is the first frame on which the
% flow vectors would be plotted. use an empty matrix '[]' for no image.
% rSize is the size of the region in which one vector is visible. scale
% over-rules the auto scaling.
imagesc(imgOriginal);colormap gray;hold on;
scale=.5;
%rSize=5;
% % Enhance the quiver plot visually by showing one vector per region
% for i=1:size(u,1)
%     for j=1:size(u,2)
%         if floor(i/rSize)~=i/rSize || floor(j/rSize)~=j/rSize
%             u(i,j)=0;
%             v(i,j)=0;
%         end
%     end
% end
%u(1,1)= -scale/1;v(1,1) = u(1,1);
u(end,1) = max(sqrt(u(:).^2+v(:).^2));v(end,1) = 0;
%u(end-1) = -u(end);v(end-1) = -v(end);
quiver(u,v, scale, 'color', 'b', 'linewidth', 2);
set(gca,'YDir','reverse');hold off;axis image;drawnow;