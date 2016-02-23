function [traj,act, score, rawValues,VValues] = spinnerOptimumP(fname,wavenums)
%use dynamic programming to determine each state's value and optimum action
%use probabilistic transitions - 1) where each action will lead on the next
%timestep (gaussian uncertainty in time-step); 2) which action an agent chooses (requires temperature param.)
%also differs from spinnerOptimum4 because it looks only n seconds into the
%future.

dat = importdata(fname);
dat = dat(ismember(dat(:,3), wavenums),:);
for i = 1:numel(wavenums)
    dat(dat(:,3) == wavenums(i),3) = i;
end
numStates = numel(unique(dat(:,2)));
wavelen = 30;
wavegap = 10;
growtime = 5;
poppabletime = .45;
slidertime = .4;%2/numStates;
popTime = .0;
dat(:,1) = dat(:,1) + (dat(:,3)-1)*(wavelen+wavegap);
dt = .01;
totallen = (wavelen + wavegap)*numel(wavenums);

travelcost = .0; %cost of moving per vertex, 0 appears to be OK!!
movecost = .01; %cost of moving at all -- prevents agent from taking a break while moving to a far target
popcost = movecost; %cost of pop action -- discourages gratuitous pressing
spawnAccess = [poppabletime 1];

rawValues = zeros(numStates,totallen/dt);
rampVal = linspace(spawnAccess(1),spawnAccess(2),diff(spawnAccess)*growtime/dt+1);
%rampTime = linspace(spawnAccess(1)*growtime,spawnAccess(2)*growtime,diff(spawnAccess)*growtime/dt+1);
theBubble = rawValues;

for i = 1:size(dat,1)
    inds = round(dat(i,1)/dt + (1:numel(rampVal)) + poppabletime*growtime/dt);
    rawValues(dat(i,2)+1,inds) = rampVal;
    theBubble(dat(i,2)+1,inds) = i;
end
VValues = zeros(numStates,size(rawValues,2)+1);
QValues = zeros(numStates,numStates*2,size(rawValues,2)+1);
isPopped = cell(size(VValues));
act = zeros(size(rawValues));


uncertainT = .1; %std / sec;
thresh = .001;
xRange = sqrt(-log(thresh)*uncertainT*growtime);
xs = (round(-xRange/dt):round(xRange/dt))*dt;
xs = xs(2:end); %THIS NEEDS TO BE DONE FOR CORRECT FILTERING
timeFilt = zeros(numel(xs),growtime/dt);
%weber's law used to estimate temporal uncertainty of different moments
%into the future. std scales linearly with time
for i = 1:size(timeFilt,2)
    timeFilt(:,i) = exp(-(xs.^2)/(uncertainT*i*dt).^2/2);
    timeFilt(:,i) = timeFilt(:,i)/sum(timeFilt(:,i));
end
for i = 1:(totallen-wavegap)/dt
    %take chunk of stimuli, blurred with temporal uncertainty
    stimPresent = rawValues(:,i) > 0;
    stimChunk = rawValues(:,i+(1:(growtime/dt+numel(xs)/2))-1);
    for j = 1:size(rawValues,1)
        if ~stimPresent(j)
            stimChunk(j,:) = 0;
        else
            ind = find(stimChunk(j,:) == 0,1);
            stimChunk(j,ind:end) = 0;
            temp = toeplitz(zeros(1,numel(xs)),stimChunk(j,:));
            stimChunk(j,1:growtime/dt) = sum(temp(:,numel(xs)/2+1:end).*timeFilt);
        end
    end
    Vlocal = zeros(numStates,growtime/dt+1);
    for ii = growtime/dt:-1:1
        for j = 1:numStates %%THIS IS WHERE I AM
            for k = 1:numStates
                indsX(k) = mod(j-k-1,numStates)+1;
                indsY(k) = round(ii+1+slidertime/dt*min(k,numStates-k));
                indsY(k+numStates) = round(indsY(k)+popTime/dt);
                indsY([k k+numStates]) = min(size(Vlocal,2),indsY([k k+numStates]));
                valNext(k) = Vlocal(indsX(k),indsY(k))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates);
                valNext(k+numStates) = Vlocal(indsX(k),indsY(k+numStates))-travelcost*min(k,numStates-k) ...
                    - movecost*(k ~= numStates)-popcost + stimChunk(j,ii);
                if ~ismember(theBubble(j,ii),isPopped{indsX(k),indsY(k+numStates)})
                    valNext(k+numStates) = valNext(k+numStates) - ;
                end
            end
            QValues(j,:,i) = valNext;
            [VValues(j,i),act(j,i)] = max(valNext);
            if act(j,i) > numStates
                isPopped{j,i} = union(theBubble(j,i),isPopped{indsX(act(j,i)-numStates),indsY(act(j,i))});
            else
                isPopped{j,i} = isPopped{indsX(act(j,i)),indsY(act(j,i))};
            end
        end
    end
end
return
for i = size(rawValues,2):-1:1
    for j = 1:numStates
        for k = 1:numStates
            indsX(k) = mod(j-k-1,numStates)+1;
            indsY(k) = round(i+1+slidertime/dt*min(k,numStates-k));
            indsY(k+numStates) = round(indsY(k)+popTime/dt);
            indsY([k k+numStates]) = min(size(VValues,2),indsY([k k+numStates]));
            valNext(k) = VValues(indsX(k),indsY(k))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates);
            valNext(k+numStates) = VValues(indsX(k),indsY(k+numStates))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates)-popcost;
            if ~ismember(theBubble(j,i),isPopped{indsX(k),indsY(k+numStates)})
                valNext(k+numStates) = valNext(k+numStates) + rawValues(j,i);
            end
        end
        QValues(j,:,i) = valNext;
        [VValues(j,i),act(j,i)] = max(valNext);
        if act(j,i) > numStates
            isPopped{j,i} = union(theBubble(j,i),isPopped{indsX(act(j,i)-numStates),indsY(act(j,i))});
        else
            isPopped{j,i} = isPopped{indsX(act(j,i)),indsY(act(j,i))};
        end
    end
end
traj = zeros(1,size(rawValues,2));
acts = traj;
traj(1) = 1;
score = zeros(size(traj));
i = 2;
while i <= size(traj,2)
    acta = act(traj(i-1),i-1);acts(i) = acta;
    score(i) = score(i-1);
    if acta > numStates
        score(i) = score(i) + rawValues(traj(i-1),i-1);
    end
    if acta > numStates
        acta = acta - numStates;
        iOld = i;
        i = i + round(popTime/dt);
        score(iOld:i) = score(iOld);
        traj(iOld:i) = traj(iOld-1);
    end
    if acta == numStates
        traj(i) = traj(i-1);%
    else
        iOld = i;
        slideTime = round(slidertime/dt*min(acta,numStates-acta));
        i = i + slideTime;
        score(iOld:i) = score(iOld)-movecost - travelcost*min(acta,numStates-acta);
        if acta < numStates-acta
            traj(iOld:i) = mod(linspace(traj(iOld-1),traj(iOld-1)-acta,slideTime+1)-.5,numStates)+.5;
        else
            traj(iOld:i) = mod(linspace(traj(iOld-1),traj(iOld-1)+numStates-acta,slideTime+1)-.5,numStates)+.5;
        end
    end
    i = i+1;
end
figure;subplot(311);%imagesc(rawValues);hold all;plot(traj,'r','linewidth',2);scatter(find(acts > numStates),traj(acts > numStates),'r','filled');
plotcirc(traj(2:end),rawValues,acts(2:end)>numStates,'k');
subplot(312);imagesc(bsxfun(@minus,VValues,min(VValues)));
[~,y] = meshgrid(1:size(act,2),1:size(act,1));
subplot(313);imagesc(mod(y-act-1,numStates)+1);
%figure;hist(act(:),0:numStates*2);
figure;plot(VValues');hold all;plot(score);
end

