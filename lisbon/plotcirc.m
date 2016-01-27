function plotcirc(traj,im,acts)
cols = ['b','r','g'];
if exist('im','var')
    imagesc(im);
end
hold all;
for i = 1:size(traj,1)
    brk = find(abs(diff(traj(i,:))) > 1);% | diff(pos(:,9) < 0));
    brk = [0 brk size(traj,2)];
    numel(brk)
    %posa = traj(i,:);%pos(:,9)  + (pos(:,3)-1)*40;%pos(:,1);
    for m = 1:numel(brk)-1
        plot((brk(m)+1):brk(m+1),traj(i,(brk(m)+1):brk(m+1)),cols(i),'linewidth',2);%posa((brk(m)+1):brk(m+1)),
    end
    %temp = dataOut(:,9)+(dataOut(:,3)-1)*40;%dataOut(:,1);
    drawnow;
    %figure(2);hold all;plot(temp,dataOut(:,4)+dataOut(:,5),cols(find(eq)));drawnow;
end
