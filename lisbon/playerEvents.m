function data = playerEvents(data,fname)
%% gets the ordered sequence of target interventions for each user
%% adds this to 'gameEvents' field in data
%gameEvents:
% 1 - event (bank/cash/burst
% 2 - wave
% 3 - score gain/loss
% 4 - which spawnEvent corresponding
% 5 - slider location
% 6 - lives lost so far
% 7 - discrepancy spawnTime / gameTime estimate
% 8 - event time within wave
% 9 - score as of last play

spawnDat = importdata(fname);
numStates = numel(unique(spawnDat(:,2)));
growTime = 5;
maxScore = 10;
anglePad = .4;
keys = {37,38,39,40}; %keyboard codes corresponding to [<,^,>,v];
%keys = {[14 18],[1 16], 15, 0}; %gamepad codes corresponding to the same (left = 14/18, bank = 1/16)

user = 1:numel(data.users);
gameVersion = zeros(1,numel(data.users));
gameVersion(2:3) = 1;
empDiff = [0 0 0 0]; %manual correction


for i = user
    session = 1:numel(data.users{i}.session);
    for j = session
        f = find(diff([0; sum(data.users{i}.session{j}.event(:,4:5),2)]) > 0);
        if gameVersion(i) == 1
            fa = f - 1;
            fa = sort([fa; find(data.users{i}.session{j}.event(:,2) == .3)]);
            f = sort([f; find(data.users{i}.session{j}.event(:,2) == .3)]);
            gameEvents = data.users{i}.session{j}.event(fa,2:end);
            gameEvents(:,3:4) = data.users{i}.session{j}.event(f,4:5);
        else
            f = sort([f; find(data.users{i}.session{j}.event(:,2) == .3)]);
            gameEvents = data.users{i}.session{j}.event(f,2:end);
        end
        gameEvents(ismember(gameEvents(:,1),keys{4}),1) = 1;
        gameEvents(ismember(gameEvents(:,1),keys{2}),1) = 2;
        gameEvents(gameEvents(:,1) == .3,1) = -1;
        gameEvents(:,6) = diff([0;gameEvents(:,8)]) < 0 & diff([0; gameEvents(:,2)]) == 0;
        gameEvents(:,6) = cumsum(gameEvents(:,6));
        temp = gameEvents(:,3);
        for k = 0:max(gameEvents(:,6))
            fRound = find(gameEvents(:,6) == k);
            gameEvents(fRound,3) = diff([temp(max(1,fRound(1)-1));sum(gameEvents(fRound,3:4),2)]);
        end
        gameEvents(:,3) = gameEvents(:,3)./abs(gameEvents(:,1));
        gameEvents(:,5) = mod(180-gameEvents(:,5),360)/360*numStates;
        if sum(~ismember(gameEvents(:,1),[-1 1 2]) > 0)
            [i j]
            return
        end
        for k = 1:size(gameEvents,1)
            if gameEvents(k,1) == -1
                estSpawn = gameEvents(k,8) - growTime - empDiff(i);
                temp = find(gameEvents(k,2) == spawnDat(:,3));
            else
                estSpawn = gameEvents(k,8) - gameEvents(k,3)/maxScore*growTime - empDiff(i);
                angDiff = gameEvents(k,5) - spawnDat(:,2);
                angDiff = min(mod(angDiff,numStates),mod(-angDiff,numStates)) < anglePad;
                temp = find(gameEvents(k,2) == spawnDat(:,3) & angDiff);
            end
            [m,in] = min(abs(estSpawn-spawnDat(temp)));
            gameEvents(k,4) = temp(in);%find(gameEvents(k,2) == spawnDat(:,3)); & abs(gameEvents(k,2)-spawnDat(:,1))<timePad);
            gameEvents(k,7) = estSpawn-spawnDat(temp(in));
        end
        %subplot(4,8,8*(i-1)+j);hist(mod(gameEvents(:,5),1),0:.01:1);
        %hist(hist(gameEvents(:,4),1:max(gameEvents(:,4))),0:2);%plot(gameEvents(:,7));%plot(diff(gameEvents(:,4)));%
        data.users{i}.session{j}.gameEvents = gameEvents;
    end
end
for i = user
%    subplot(2,2,i);
end


% for i = 1:numStates
%     temp = dat(dat(:,3) == i,:);
%     minDist = min(mod(diff(temp(:,2)),numStates),mod(-diff(temp(:,2)),numStates));
%     subplot(2,3,i);
%     %% 2-body
%     %scatter(minDist,diff(temp(:,1)));hold all;
%     %% 2-body (correcting for travel time)
%     %scatter(minDist,diff(temp(:,1))-minDist*sliderTime);
%     %% 3-body
%     
%     %% 3-body (correcting for travel time)
% end
