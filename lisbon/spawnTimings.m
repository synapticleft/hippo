function [total_objs,event_log] = spawnTimings(fname,n_positions)
%% generates a text file with wave number, location, and appearance time
%% currently, rate of appearance is increased in each wave, and bubbles are
%% equally likely to appear at each site. Wave lengths are currently all 10 s.
%% inputs - fname is the name of the file to be created, n_positions is the # of vertices (4, 6, or 8)
%% outputs the text file in current directory

%traversal_time = .5; %in seconds, for slider
%gap_between_waves = 5; %in seconds
wave_length = 30; %in seconds
n_waves = 10; %number of waves
%n_positions = 8; %number of places object can appear
appearance_rate = linspace(1,5,n_waves)/6;%n_positions; %per second
grow_time = 6; %in seconds, time for which slider overlaps with spawn
dt = .01;%traversal_time; %discretization of time
event_log = [];
counter = 1;
total_objs = zeros(n_waves,2);

for i = 1:n_waves
    events = min(1,poissrnd(appearance_rate(i)*dt,n_positions,wave_length/dt));
    total_objs(i,1) = sum(events(:));
     for j = 1:size(events,2)
         for k = 1:size(events,1)
             if events(k,j)
                 inds = j+(1:grow_time/dt);
                 inds(inds > size(events,2)) = []; 
                 %events(k,inds) = 0; %makes sure another bubble doesnt appear in the future while current one still there
                 event_log(counter,:) = [j*dt k-1 i];% + (wave_length+gap_between_waves)*(i-1) %if we want time to accumulate across waves
                 counter = counter + 1;
                 total_objs(i,2) = total_objs(i,2) + 1;
             end
         end
     end
end
plot(event_log);
%dlmwrite(fname,event_log);

%total_objs % outputs total events in each wave, before and after correcting for multiple appearances at the same spot.