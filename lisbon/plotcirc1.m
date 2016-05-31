function plotcirc1(traj,im,acts,cols)
% allow you to plot lines of the player's trajectory, but since the track
% is circular, interrupts the plot when player goes from vertex '6' to '1' 
% to prevent crowding. plots each player in a different color, as well as
% their actions as triangles.

% INPUTS
% traj = the output of getKeyStrokesConstruct
% im = matrix of the times and locations of the bubbles, serves as
%   a background
% acts = used only by plotcirc (legacy format)
% cols = colors you want to assign to each player (has a default setting)

% If we want to compare the internal consistency of timing as measured by
% time reported by the game, to the time as estimated by the score upon a
% bank event, set this to 1 to generate some plots
scoreVsTimeCheck = 0; 

% timestep, needs to match that used to generate 'im' (done by spinnerOptimum)
dt = .05;

%number of vertices, needs to match the trajectories
nVertices = 6;

if scoreVsTimeCheck == 1
f1 = figure;
f2 = figure;
offShifts = [];
end

%keys = {37,38,39,40}; %keyboard codes corresponding to [<,^,>,v];
keys = {[14],[1], 15, 0}; %gamepad codes corresponding to the same (left = 14/18, bank = 1/16)
shiftrange = -1:dt:1;
if ~exist('cols','var')
    cols = ['b','r','g','y'];
end
if exist('im','var') && ~isempty(im)
    imagesc((1:size(im,2))*dt,(1:size(im,1)),im);
end
hold all;
if isstruct(traj) %getKeyStrokesConstruct used to generate data structure containing game events
    symbols = ['<','^','>','v'];
    for j = 1:numel(traj.users)
        for i = 1:numel(traj.users{j}.session)
            events = traj.users{j}.session{i}.gameData.event;
            events(events(:,3) ~= 5,:) = [];
            events(:,3) = 1;
            time = events(:,9)+(events(:,3)-1)*40;%fullData.users{j}.session{i}.position(:,1);
            position = mod(180-events(:,6),360)/360*nVertices + 1;% + 1;
            if scoreVsTimeCheck == 1
                inds = find(ismember(events(:,2),keys{2}));
                values = max(0,events(inds,5) - events(inds-1,5));
                offset = zeros(numel(inds),numel(shiftrange));
                for n= 1:numel(inds)
                    offset(n,:) = abs(im(mod(round(position(inds(n)))-1,nVertices)+1,round((time(inds(n))-shiftrange)/dt))-values(n)/20);
                end
                offset(offset > .4 | offset == 0) = nan;
                [~,offShift] = min(offset');
                offShifts = [offShifts offShift];
                figure(f1);hold all;plot(shiftrange,offset,cols(j));
                figure(f2);
            end
            hold all;
            for m = 1:4
                inds = ismember(events(:,2),keys{m});% == keys(m);
                scatter(time(inds),position(inds),cols(j),symbols(m),'filled');%-(offShift-1)*dt-min(shiftrange)
            end
            brk = find(abs(diff(position)) > 1 | diff(events(:,9)) < 0);
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