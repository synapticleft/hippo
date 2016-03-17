function [event_log] = spawnTimings1(fname,n_positions,slide_time)
%% generates a text file with wave number, location, and appearance time
%% In this version, each wave is comprised of subwaves, each of which consists of 
%% between 2 and n_positions targets, spawned within a fixed window of time.
%% Larger waves contain more targets per subwave.
%% outputs the text file in current directory

n_subwaves = 20;
n_waves = n_positions - 1;
%spawn_window = 2; %time within which all spawns occur in a subwave
spawn_window = slide_time*maxSpawnTime(n_positions);
gap_subwaves = 3; %gap between subwaves
dt = .01; %discretization of time

% wave_length = 30; %in seconds
% n_waves = 10; %number of waves
% %n_positions = 8; %number of places object can appear
% appearance_rate = linspace(1,5,n_waves)/6;%n_positions; %per second
% grow_time = 6; %in seconds, time for which slider overlaps with spawn THIS IS TOO HIGH!!!!!!!!!!!!!!!!!!!!!!
% event_log = [];
% counter = 1;
% total_objs = zeros(n_waves,2);

counter = 0;
for i = 1:n_waves
    for j = 1:n_subwaves
        r = randn([round(spawn_window(i)/dt) i]);
        order = randperm(n_positions);
        event_log(counter+1,1) = 1; %first spawn is at start of subwave
        [~,event_log(counter+(2:i+1),1)] = max(r);
        event_log(counter+(1:i+1),1) = sort(event_log(counter+(1:i+1),1))*dt + (j-1)*(gap_subwaves+spawn_window(i));
        event_log(counter+(1:i+1),2) = order(1:i+1)-1;
        event_log(counter+(1:i+1),3) = i;
        counter = counter + i+1;
    end
end
dlmwrite(fname,event_log);

function t = maxSpawnTime(nVertices)
t = 0;
for i = 1:nVertices - 1
    if mod(i,2)
        t(i+1) = t(i) + ceil((nVertices-1)/2);
    else
        t(i+1) = t(i) + floor((nVertices-1)/2);
    end
end
t = t(2:end);

%     events = min(1,poissrnd(appearance_rate(i)*dt,n_positions,wave_length/dt));
%     total_objs(i,1) = sum(events(:));
%      for j = 1:size(events,2)
%          for k = 1:size(events,1)
%              if events(k,j)
%                  inds = j+(1:grow_time/dt);
%                  inds(inds > size(events,2)) = []; 
%                  %events(k,inds) = 0; %makes sure another bubble doesnt appear in the future while current one still there
%                  event_log(counter,:) = [j*dt k-1 i];% + (wave_length+gap_between_waves)*(i-1) %if we want time to accumulate across waves
%                  counter = counter + 1;
%                  total_objs(i,2) = total_objs(i,2) + 1;
%              end
%          end
%      end
% end
% plot(event_log);


%total_objs % outputs total events in each wave, before and after correcting for multiple appearances at the same spot.