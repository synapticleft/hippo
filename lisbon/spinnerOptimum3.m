function [traj, score, rawValues] = spinnerOptimum3(fname,wavenums)
%use dynamic programming to determine each state's value and optimum action

%r =  [.1 .2 .3 .4 .5 .4 .3 .2];%[.1 .3 .1 .2 .1 .4 .5 .6];

dat = importdata(fname);
dat = dat(ismember(dat(:,3), wavenums),:);
wavelen = 30;
wavegap = 10;
growtime = 5;
poptime = .45;
slidertime = .1;
dat(:,1) = dat(:,1) + (dat(:,3)-1)*(wavelen+wavegap);
numStates = numel(unique(dat(:,2)));
dt = .01;
totallen = (wavelen + wavegap)*numel(wavenums);

travelcost = .01;
spawnAccess = [poptime 1];

rawValues = zeros(numStates,totallen/dt);
rampVal = linspace(spawnAccess(1),spawnAccess(2),diff(spawnAccess)*growtime/dt+1);
theBubble = rawValues;

for i = 1:size(dat,1)
    inds = round(dat(i,1)/dt + (1:numel(rampVal)));
    rawValues(dat(i,2)+1,inds) = rampVal;
    theBubble(dat(i,2)+1,inds) = i;
end
dpValues = zeros(size(rawValues,1),size(rawValues,2)+1);
isPopped = cell(size(dpValues));
act = zeros(size(rawValues));
for i = size(rawValues,2):-1:1
    for j = 1:numStates
        for k = 1:numStates
            indsX(k) = mod(j-k-1,numStates)+1;
            indsY(k) = round(min(size(dpValues,2),i+1+slidertime/dt*min(k,numStates-k)));
%        end
        %indsX = [mod(j-1-1,numStates)+1  j mod(j+1-1,numStates)+1];
        %indsY = [round(min(size(dpValues,2),i+1+slidertime/dt))  i+1 round(min(size(dpValues,2),i+1+slidertime/dt))];
%        for k = 1:3
            valNext(k) = dpValues(indsX(k),indsY(k))-travelcost*min(k,numStates-k);
            %if k ~= numStates
            %    valNext(k) = valNext(k) - travelcost;
            %end
            if ~ismember(theBubble(j,i),isPopped{indsX(k),indsY(k)})
                valNext(k+numStates) = valNext(k) + rawValues(j,i);
            else
                valNext(k+numStates) = -Inf;
            end
        end
        [dpValues(j,i),act(j,i)] = max(valNext);
        if act(j,i) > numStates
            isPopped{j,i} = union(theBubble(j,i),isPopped{indsX(act(j,i)-numStates),indsY(act(j,i)-numStates)});
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
        acta = acta - numStates;
    end
    if acta == numStates
        traj(i) = traj(i-1);%
    else
        iOld = i;
        i = i + slidertime/dt*min(acta,numStates-acta);
        score(iOld:i) = score(iOld);
        if acta < numStates-acta
            traj(iOld:i) = mod(linspace(traj(iOld-1),traj(iOld-1)-acta,slidertime/dt*min(acta,numStates-acta)+1)-.5,numStates)+.5;
        else
            traj(iOld:i) = mod(linspace(traj(iOld-1),traj(iOld-1)+numStates-acta,slidertime/dt*min(acta,numStates-acta)+1)-.5,numStates)+.5;
        end
    end
    i = i+1;
end
figure;plot(score);
figure;plot(traj);
figure;subplot(311);imagesc(rawValues);hold all;plot(traj,'r','linewidth',2);scatter(find(acts > numStates),traj(acts > numStates),'r','filled');
subplot(312);imagesc(dpValues);
subplot(313);imagesc(act);
figure;hist(act(:),0:numStates*2);
figure;plot(dpValues');
end

