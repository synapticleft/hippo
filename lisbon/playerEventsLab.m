function data = playerEventsLab(data,fname)
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
%% to do: is there a way of inferring which targets were intended but missed?
%% occasionally a bubble vanishes at the very end of wave! LET TIAGO KNOW

spawnDat = importdata(fname);
numStates = numel(unique(spawnDat(:,2)));
growTime = 5;
maxScore = 10;
anglePad = .49;
%keys = {37,38,39,40}; %keyboard codes corresponding to [<,^,>,v];
%keys = {[14 18],[1 16], 15, 0}; %gamepad codes corresponding to the same (left = 14/18, bank = 1/16)
keys = {14,1,15,0};

%gameVersion = [0 0 0];%zeros(1,numel(data.users));
%gameVersion(2:3) = 1;
%empDiff = [0 0 0 0]; %manual correction

for i = 1:numel(data.users)
    for j = 1:numel(data.users{i}.session)
%        f = find(diff([0; sum(data.users{i}.session{j}.mergeData.merge_events(:,4:5),2)]) > 0);
%        if gameVersion(i) == 1
 %           fa = f - 1;
 %           fa = sort([fa; find(data.users{i}.session{j}.event(:,2) == .3)]);
 %           f = sort([f; find(data.users{i}.session{j}.event(:,2) == .3)]);
 %           gameEvents = data.users{i}.session{j}.event(fa,2:end);
 %           gameEvents(:,3:4) = data.users{i}.session{j}.event(f,4:5);
%        else
%            f = sort([f; find(data.users{i}.session{j}.event(:,2) == .3)]);
%            gameEvents = data.users{i}.session{j}.event(f,2:end);
%        end
        mergeEvents = data.users{i}.session{j}.mergeData.merge_events(:,2:end);
        changeScore = diff([0;sum(mergeEvents(:,3:4),2)]);
        removeInds = mergeEvents(:,1) == .1 | mergeEvents(:,1) == .2 | ismember(mergeEvents(:,1),keys{1}) | ismember(mergeEvents(:,1),keys{3});
        mergeEvents(removeInds,:) = [];
        changeScore(removeInds) = [];
        gameEvents = zeros(size(mergeEvents,1),10);
        gameEvents(ismember(mergeEvents(:,1),keys{4}),1) = 1;
        gameEvents(ismember(mergeEvents(:,1),keys{2}),1) = 2;
        %gameEvents(ismember(mergeEvents(:,1),keys{1}),1) = 3;
        %gameEvents(ismember(mergeEvents(:,1),keys{3}),1) = 4;
        removeInds = changeScore == 0 & gameEvents(:,1) > 0;
        gameEvents(removeInds,:) = [];
        mergeEvents(removeInds,:) = [];
        gameEvents(mergeEvents(:,1) == .3,1) = -1;
        gameEvents(:,[2 8]) = mergeEvents(:,[2 6]);
        gameEvents(:,6) = diff([0;mergeEvents(:,6)]) < 0 & diff([0; mergeEvents(:,2)]) == 0;
        gameEvents(:,6) = cumsum(gameEvents(:,6));
%        temp = gameEvents(:,3);
        for k = 0:max(gameEvents(:,6))
            fRound = find(gameEvents(:,6) == k);
            gameEvents(fRound,3) = diff([mergeEvents(max(1,fRound(1)-1),3);sum(mergeEvents(fRound,3:4),2)]);
        end
        gameEvents(:,3) = gameEvents(:,3)./abs(gameEvents(:,1));
        gameEvents(:,5) = mod(180-mergeEvents(:,5),360)/360*numStates;
        for k = 1:size(gameEvents,1)
            if gameEvents(k,1) > 0
                estSpawn = gameEvents(k,8) - gameEvents(k,3)/maxScore*growTime;% - empDiff(i);
                angDiff = gameEvents(k,5) - spawnDat(:,2);
                angDiff = min(mod(angDiff,numStates),mod(-angDiff,numStates)) < anglePad;
                temp = find(gameEvents(k,2) == spawnDat(:,3) & angDiff);
            end
            [m,in] = min(abs(estSpawn-spawnDat(temp,1)));
            gameEvents(k,4) = temp(in);%find(gameEvents(k,2) == spawnDat(:,3)); & abs(gameEvents(k,2)-spawnDat(:,1))<timePad);
            gameEvents(k,7) = estSpawn-spawnDat(temp(in),1);
        end
        f = find(gameEvents(:,1) == -1);
        for k = 1:numel(f)%find(gameEvents(:,1) == -1)
                estSpawn = gameEvents(f(k),8) - growTime;% - empDiff(i);
                temp = find(gameEvents(f(k),2) == spawnDat(:,3) & ~ismember(1:size(spawnDat,1),gameEvents(:,4))');
                [m,in] = min(abs(estSpawn-spawnDat(temp,1)));
                gameEvents(f(k),4) = temp(in);
 %               gameEvents(f(k),7) = estSpawn-spawnDat(temp(in),1);
        end
        gameEvents(f,4) = sort(gameEvents(f,4));
        gameEvents(f,7) = gameEvents(f,8)-growTime-spawnDat(gameEvents(f,4),1);
        %subplot(4,8,8*(i-1)+j);hist(mod(gameEvents(:,5),1),0:.01:1);
        %hist(hist(gameEvents(:,4),1:max(gameEvents(:,4))),0:2);%plot(gameEvents(:,7));%plot(diff(gameEvents(:,4)));%
        data.users{i}.session{j}.gameEvents = gameEvents;
    end
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
