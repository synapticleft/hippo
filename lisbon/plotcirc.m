function plotcirc(traj,im,acts,cols)

% f1 = figure;
% f2 = figure;
offShifts = [];
dt = .01;
nVertices = 6;
shiftrange = -1:dt:1;
if ~exist('cols','var')
    cols = ['b','r','g','y'];
end
if exist('im','var') && ~isempty(im)
    imagesc((1:size(im,2))*dt,(1:size(im,1)),im);
end
hold all;
if isstruct(traj)
    symbols = ['<','^','>','v'];
    for j = 1:numel(traj.users)
        for i = 1:numel(traj.users{j}.session)
            time = traj.users{j}.session{i}.event(:,9)+(traj.users{j}.session{i}.event(:,3)-1)*40;%fullData.users{j}.session{i}.position(:,1);
            position = mod(180-traj.users{j}.session{i}.event(:,6),360)/360*nVertices + 1;% + 1; 
                inds = find(traj.users{j}.session{i}.event(:,2) == 38);
                values = max(0,traj.users{j}.session{i}.event(inds,5) - traj.users{j}.session{i}.event(inds-1,5));
                offset = zeros(numel(inds),numel(shiftrange));
                for n= 1:numel(inds)
                    offset(n,:) = abs(im(mod(round(position(inds(n)))-1,nVertices)+1,round((time(inds(n))-shiftrange)/dt))-values(n)/20);
                end
                offset(offset > .4 | offset == 0) = nan;
                [~,offShift] = min(offset');
                offShifts = [offShifts offShift];
%                  figure(f1);hold all;plot(shiftrange,offset,cols(j));
%              figure(f2);hold all;
            for m = [2 4]
                inds = traj.users{j}.session{i}.event(:,2) == m+36;
                scatter(time(inds),position(inds),cols(j),symbols(m),'filled');%-(offShift-1)*dt-min(shiftrange)
            end
            brk = find(abs(diff(position)) > 1 | diff(traj.users{j}.session{i}.event(:,9)) < 0);
            brk = [0 brk'];
            for m = 1:numel(brk)-1
                plot(time((brk(m)+1):brk(m+1)),position((brk(m)+1):brk(m+1)),cols(j),'linewidth',2);%posa((brk(m)+1):brk(m+1)),-(offShift-1)*dt-min(shiftrange)
            end
        end
    end
else
    for i = 1:size(traj,1)
        brk = find(abs(diff(traj(i,:))) > 1);% | diff(pos(:,9) < 0));
        brk = [0 brk size(traj,2)];
        %posa = traj(i,:);%pos(:,9)  + (pos(:,3)-1)*40;%pos(:,1);
        for m = 1:numel(brk)-1
            plot(((brk(m)+1):brk(m+1))*dt,traj(i,(brk(m)+1):brk(m+1)),'k','linewidth',2);%posa((brk(m)+1):brk(m+1)),
            
        end
        if exist('acts','var') && ~isempty(acts)
            scatter(find(acts(i,:))*dt,traj(i,acts(i,:)),cols,'filled');
        end
        %temp = dataOut(:,9)+(dataOut(:,3)-1)*40;%dataOut(:,1);
        %drawnow;
        %figure(2);hold all;plot(temp,dataOut(:,4)+dataOut(:,5),cols(find(eq)));drawnow;
    end
end
% figure;plot(offShifts)
% figure;hist(shiftrange(offShifts),-.05:.01:1);