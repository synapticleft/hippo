function [VValues,QValues,PValues,traj,acts, score] = spinnerOptimumPData(fname,VValues,QValues)
%use dynamic programming to determine each state's value and optimum action
%use probabilistic transitions - 1) where each action will lead on the next
%timestep (gaussian uncertainty in time-step); 2) which action an agent chooses (requires temperature param.)
%also differs from spinnerOptimum4 because it looks only n seconds into the
%future.

temperature = .0001;%.01;%temperature for soft-max action selection (larger than .0015 explodes to inf)
penalty = 2; % amount that score decreases due to ability to pop target in the future.

growtime = 5;
poppabletime = .45;
dt = .05;
slidertime = 2/6;%dt;%2/numStates;
popTime = .0;
allDat = importdata(fname);
numStates = numel(unique(allDat(:,2)));

timecost = 0*dt; %cost of moving per vertex, 0 appears to be OK!!
movecost = 0.01; %cost of moving at all -- prevents agent from taking a break while moving to a far target
popcost = 0; %cost of pop action -- discourages gratuitous pressing
spawnAccess = [poppabletime 1];
rampVal = linspace(spawnAccess(1),spawnAccess(2),diff(spawnAccess)*growtime/dt);

for l = 1:5
dat = allDat(ismember(dat(:,3),l),:);
dat(:,3) = 1;
dat(:,1) = round(dat(:,1)/dt)*dt;
totallen = 35;

rawValues = zeros(numStates,totallen/dt);

for i = 1:size(dat,1)
    inds = round((dat(i,1) + growtime)/dt);
    rawValues(dat(i,2)+1,(inds-numel(rampVal)+1):inds) = rampVal;
end
VValues = zeros(numStates,size(rawValues,2));
QValues = zeros(numStates,numStates*2,size(rawValues,2));
PDensity = zeros(size(VValues));
PDensity(:,1) = 1/numStates;
%isPopped = cell(size(VValues));
%act = zeros(size(rawValues));

uncertainT = .01; %std / sec;
thresh = .01;
xRange = sqrt(-log(thresh)*uncertainT*growtime);
xs = (round(-xRange/dt):round(xRange/dt))*dt;
xPad = (numel(xs)-1)/2;
%xs = xs(2:end); %THIS NEEDS TO BE DONE FOR CORRECT FILTERING
timeFilt = zeros(numel(xs),growtime/dt+xPad);
%weber's law used to estimate temporal uncertainty of different moments
%into the future. std scales linearly with time
for i = 1:size(timeFilt,2)
    timeFilt(:,i) = exp(-(xs.^2)/(uncertainT*i*dt).^2/2);
    timeFilt(:,i) = timeFilt(:,i)/sum(timeFilt(:,i));
end
%stimChunk = zeros(numStates,growtime/dt+numel(xs)/2);
rawChunk = zeros(1,growtime/dt+numel(xs)-1);%/2
Vlocal = zeros(numStates,growtime/dt+xPad+1);
Qlocal = zeros(numStates,numStates*2,size(Vlocal,2));
valNext = zeros(1,2*numStates);
inds = [];
actTime = repmat(round(max(1,slidertime/dt*min((1:numStates),numStates-(1:numStates)))),[1 2]);
%temperatures = temperature*actTime;
for i = 1:size(rawValues,2)
    %take chunk of stimuli, blurred with temporal uncertainty
    rawChunk(:) = 0;
    visBubbles = find(round((dat(:,1)+growtime)/dt) >= i & round(dat(:,1)/dt) <= i);%%dat(:,1) > i*dt - growtime & dat(:,1) <= i*dt);
    if numel(visBubbles)
        stimChunk = zeros(numel(visBubbles),growtime/dt+xPad);
        for k = 1:numel(visBubbles) %  
            inds(2) = round((dat(visBubbles(k),1)+growtime)/dt - i+1);%DOESNT WORK FOR WEIRD DT.. FIX; ALSO i MINUS 1
            inds(1) = max(1,inds(2)-numel(rampVal)+1);%
%            inds = max(round((dat(visBubbles(k),1) + growtime*poppabletime)/dt-i+2),1):round((dat(visBubbles(k),1)+growtime)/dt - i+1); 
            rawChunk(:) = 0;
            rawChunk(inds(1):inds(2)) = rampVal(end-diff(inds):end);
            temp = toeplitz(rawChunk,zeros(1,numel(xs)))';
            stimChunk(k,:) = sum(temp(:,xPad+1:end).*timeFilt); %(:,numel(xs)/2+1:end)
            %         rawChunk(:) = 0; %%THIS MAY BE USEFUL IN ORDER TO INCORPORATE BANKING
            %         rawChunk(inds) = 1;
            %         temp = toeplitz(zeros(1,numel(xs)),rawChunk);
            %         stimChunk(k,:) = sum(temp(:,numel(xs)/2+1:end).*timeFilt);
        end
        Vlocal(:) = 0;
        Qlocal(:) = 0;
        VBlocal = zeros(numStates,numel(visBubbles),size(Vlocal,2));
        VBtemp = zeros(numStates*2,numel(visBubbles));
        for ii = size(stimChunk,2):-1:1
            for j = 1:numStates
                curBubble = find(dat(visBubbles,2)+1 == j);
                indsX = mod(j-(1:numStates)-1,numStates)+1;
                indsY = min(size(Vlocal,2),ii+actTime);
                for k = 1:numStates
                    %indsX(k) = mod(j-k-1,numStates)+1;
                    %indsY([k k+numStates]) = min(size(Vlocal,2),round(ii+max(1,slidertime/dt*min(k,numStates-k))));
                    %indsY(k+numStates) = min(size(Vlocal,2),round(indsY(k)+popTime/dt)); %IF POPTIME IS NONZERO
                    %indsY([k k+numStates]) = min(size(Vlocal,2),indsY([k k+numStates]));
                    VBtemp(k,:) = VBlocal(indsX(k),:,indsY(k)); %Needed??
                    VBtemp(k+numStates,:) = VBtemp(k,:);
                    valNext([k k+numStates]) = Vlocal(indsX(k),indsY(k))-timecost*(actTime(k)-1) - movecost*(k ~= numStates); %SINCE POPTIME = 0 NOW
                     if numel(curBubble) && stimChunk(curBubble,ii) > 0
                        VBtemp(k+numStates,curBubble) = stimChunk(curBubble,ii);
                        %valNext(k+numStates) = Vlocal(indsX(k),indsY(k+numStates))-travelcost*min(k,numStates-k)- movecost*(k ~= numStates) ...
                        valNext(k+numStates) = valNext(k+numStates) - popcost + stimChunk(curBubble,ii) - penalty*VBlocal(indsX(k),curBubble,indsY(k+numStates));
                    end
                end
                Qlocal(j,:,ii) = valNext;
                pAct = exp(valNext./temperature);%log(1+exp(valNext./temperature); 1./(1+exp((valNext-V0)/temperature));max(0,valNext)
                %pAct = max(eps,valNext); APPEARS TO BE TOO NONSELECTIVE??
                pAct = pAct./actTime;
                if isinf(sum(pAct))
                    pAct = double(valNext == max(valNext));
                end
                pAct = pAct/sum(pAct);
                Vlocal(j,ii) = pAct*valNext';
                VBlocal(j,:,ii) = pAct*VBtemp;
                if ii == 1
                    PValues(j,:,i) = pAct;
                end
%                pActs(j,:,ii) = pAct;
            end
        end
        VValues(:,i) = Vlocal(:,1);
        QValues(:,:,i) = Qlocal(:,:,1);
%        plot(Vlocal');title(i);drawnow;
    else
        QValues(:,[1:numStates-1 numStates+(1:numStates)],i) = -movecost;%ASSUMING NO TIMECOST
        %QValues(:,numStates+(1:numStates),i) = QValues(:,1:numStates,i);
        pAct = exp(QValues(:,:,i)/temperature);%log(1+exp(valNext./temperature); 1./(1+exp((valNext-V0)/temperature));max(0,valNext)
        %pAct = max(eps,valNext); APPEARS TO BE TOO NONSELECTIVE??
        pAct = bsxfun(@rdivide,pAct,actTime);
        if isinf(sum(pAct))
            pAct = double(valNext == max(valNext));
        end
        PValues(:,:,i) = bsxfun(@rdivide,pAct,sum(pAct,2));
    end
    PDensity(:,i) = PDensity(:,i)/(sum(PDensity(:,i)));
    tempP = bsxfun(@times,PValues(:,1:numStates,i) + PValues(:,(1:numStates)+numStates,i),PDensity(:,i));
    for j = 1:numStates
        indsX = mod(j-(1:numStates)-1,numStates)+1;
        indsY = min(size(VValues,2),i+actTime);
        for k = 1:numStates
            PDensity(indsX(k),indsY(k)) = PDensity(indsX(k),indsY(k)) + tempP(j,k);
        end
    end
end
PDensity = min(PDensity,1);

traj = zeros(1,size(rawValues,2));
acts = traj;%acta = 1;
traj(1) = 1;
score = zeros(size(traj));
i = 1;%iOld = 1;
while i <= size(rawValues,2)
    [~,acta] = max(PValues(traj(i),:,i));%mod(acts(iOld)-1,numStates)+1
    if i >= 2
        score(i) = score(i-1);
    else
        score(i) = 0;
    end
    acts(i) = acta;
    if acta > numStates
        score(i) = score(i) + rawValues(traj(i),i);
        acta = acta - numStates;
        %iOld = i;
        %i = i + round(popTime/dt);
        %score(iOld:i) = score(iOld);
        %traj(iOld:i) = traj(iOld-1);
    end
    if acta == numStates
        traj(i+1) = traj(i);%
        i = i+1;
    else
%         if acta <= numStates/2 %%THIS IS NOT WORKING
%             acta = 1;
%         else
%             acta = numStates-1;
%         end
        iOld = i;
        slideTime = round(slidertime/dt*min(acta,numStates-acta));
        i = i + slideTime;
        score(iOld:i) = score(iOld)-movecost - timecost*min(acta,numStates-acta);
        if acta < numStates-acta
            traj(iOld:i) = mod(linspace(traj(iOld),traj(iOld)-acta,slideTime+1)-.5,numStates)+.5;
        else
            traj(iOld:i) = mod(linspace(traj(iOld),traj(iOld)+numStates-acta,slideTime+1)-.5,numStates)+.5;
        end
    end

    %i = i+1;
end
model.rawValues{i} = rawValues;
model.trajectory = [model.trajectory; 
end
figure;subplot(311);%imagesc(rawValues);hold all;plot(traj,'r','linewidth',2);scatter(find(acts > numStates),traj(acts > numStates),'r','filled');
plotcirc(traj,rawValues,acts>numStates,'k');
%subplot(312);imagesc(VValues);%bsxfun(@minus,VValues,min(VValues))
subplot(312);imagesc(rawValues);
hold on;
for i = 1:numStates
    [~,I] = max(squeeze(PValues(i,end:-1:1,:)));
%    f = find(squeeze(QValues(i,:,:)) == max(squeeze(QValues(i,:,:))),1,'last');
    scatter(find(I <= numStates & rawValues(i,:) > 0),i*ones(1,sum(I <= numStates & rawValues(i,:) > 0)),'k','filled');%&
end
subplot(313);imagesc(PDensity);
% totalPlayers = zeros(size(VValues));
% totalPlayers(:,1) = 1/size(totalPlayers,1);
% for i = 2:size(Pvalues,3)
%     num = squeeze(totalPlayers(:,i-1).*squeeze(PValues(:,:,i));
% end
%[~,y] = meshgrid(1:size(act,2),1:size(act,1));
%subplot(313);imagesc(mod(y-act-1,numStates)+1);
%figure;hist(act(:),0:numStates*2);
%figure;plot(VValues');%hold all;plot(score);
score(end)
end

    
