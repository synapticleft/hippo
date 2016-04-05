function [event_log] = spawnTimings3(fname,n_positions,slide_time)
%% generates a text file with wave number, location, and appearance time
%% In this version, each wave is comprised of subwaves, each of which consists of 
%% between 2 and n_positions targets, spawned within a limited time.
%% Larger waves contain more targets per subwave.
%% Uses 3 algorithms:
%% 1) Spawn times are sampled independently within the subwave period
%% 2) Spawn times are mutually exclusive and discretized (default to slide_time)
%% 3) Spawn order is chosen, and timing is determined to be the time it takes 
%%      to reach the target, plus noise sampled from a gaussian (default to .1*slide_time)
%% In each case, 1st spawn happens at time 0, and subwave is terminated after last spawn
%% outputs the text file in current directory

n_subwaves = 20;
n_waves = n_positions - 1;
spawn_window = slide_time*maxSpawnTime(n_positions);
gap_subwaves = 3; %gap between subwaves
dt = .01; %discretization of time
dt2 = slide_time;
jitter = .2;

event_log{3} = [];
counter = 0;
for i = 1:n_waves
    for j = 1:n_subwaves
        r = randn([round(spawn_window(i)/dt) i]);
        order = randperm(n_positions);
        event_log{1}(counter+1,1) = 0; %first spawn is at start of subwave
        [~,event_log{1}(counter+(2:i+1),1)] = max(r);
        event_log{1}(counter+(1:i+1),1) = sort(event_log{1}(counter+(1:i+1),1))*dt;
        event_log{1}(counter+(1:i+1),2) = order(1:i+1)-1;
        event_log{1}(counter+(1:i+1),3) = i;
        r = randn([round(spawn_window(i)/dt2) 1]);
        [~,r] = sort(r,'descend');
        event_log{2}(counter+(1:i+1),1) = [0; sort(r(1:i))]*dt2;
        event_log{2}(counter+(1:i+1),2) = order(1:i+1)-1;
        event_log{2}(counter+(1:i+1),3) = i;
        event_log{3}(counter+(1:i+1),2) = order(1:i+1)-1;
        event_log{3}(counter+(1:i+1),3) = i;  
        event_log{3}(counter+1,1) = 0;
        for k = 1:i %%I WOULD LIKE TO MAKE THIS USE THE PREVIOUS POINT THAT WOULD CAUSE THE LATEST SPAWN
            %travelTime = slide_time*min(abs(event_log{3}(counter+1+k,2)-event_log{3}(counter+k,2)),...
            %    abs(event_log{3}(counter+1+k,2)-(event_log{3}(counter+k,2)+1-n_positions)));
            travelTime = slide_time*min(mod(event_log{3}(counter+1+k,2)-event_log{3}(counter+k,2),n_positions),...
                mod(event_log{3}(counter+k,2)-event_log{3}(counter+1+k,2),n_positions));
            event_log{3}(counter+1+k,1) = event_log{3}(counter+k,1) + travelTime + randn*jitter*slide_time;
        end
        event_log{3}(counter+(1:i+1),1) = event_log{3}(counter+(1:i+1),1);
        if j > 1
            for k = 1:3
                event_log{k}(counter+(1:i+1),1) = event_log{k}(counter+(1:i+1),1) + gap_subwaves + event_log{k}(counter,1);
            end
        end      
        counter = counter + i + 1;
    end
end
for i = 1:3
    dlmwrite([fname num2str(i) '.txt'],event_log{i});
end

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