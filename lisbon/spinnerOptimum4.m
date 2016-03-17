function [dpValues,rawValues,traj,act, score,travelLength] = spinnerOptimum4(fname,wavenums,slidertime)
% uses dynamic programming to determine each state's value and optimum action
% INPUT
% fname = input file of bubble spawn times (eg. 'abs6test.txt')
% wavenums = which waves to consider (can be just one, or several, e.g. 1:10)
% OUTPUT
% dpValues = the value of being at a given vertex, at a given time
% rawValues = the value of popping a bubble at a given vertex, at a given
    % time (used by plotcirc)
% traj = a time series of the trajectory used by the agent
% acts = the actions taken by the agent
% score = the cumulative score of the agent

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
if ~exist('slidertime','var')
slidertime = .1;
end
popTime = .0;
dat(:,1) = dat(:,1) + (dat(:,3)-1)*(wavelen+wavegap);
dt = .05;
totallen = (wavelen + wavegap)*numel(wavenums);

travelcost = .0; %cost of moving per vertex, 0 appears to be OK!!
movecost = .01; %cost of moving at all -- prevents agent from taking a break while moving to a far target
spawnAccess = [poppabletime 1];

rawValues = zeros(numStates,totallen/dt);
rampVal = linspace(1/growtime,spawnAccess(2),diff(spawnAccess)*growtime/dt);%spawnAccess(1)
%rampTime = linspace(spawnAccess(1)*growtime,spawnAccess(2)*growtime,diff(spawnAccess)*growtime/dt+1);
theBubble = rawValues;

for i = 1:size(dat,1)
    inds = round(dat(i,1)/dt + (1:numel(rampVal)) + poppabletime*growtime/dt);
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
            indsY(k) = round(i+1+slidertime/dt*min(k,numStates-k));
            indsY(k+numStates) = round(indsY(k)+popTime/dt);
            indsY([k k+numStates]) = min(size(dpValues,2),indsY([k k+numStates]));
            valNext(k) = dpValues(indsX(k),indsY(k))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates);
            valNext(k+numStates) = dpValues(indsX(k),indsY(k+numStates))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates);
            if ~ismember(theBubble(j,i),isPopped{indsX(k),indsY(k+numStates)})
                valNext(k+numStates) = valNext(k+numStates) + rawValues(j,i);
            else
                valNext(k+numStates) = -Inf;
            end
        end
        [dpValues(j,i),act(j,i)] = max(valNext);
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
travelLength = 0;
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
        travelLength = travelLength + min(acta,numStates-acta);
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
subplot(312);imagesc(bsxfun(@minus,dpValues,min(dpValues)));
[~,y] = meshgrid(1:size(act,2),1:size(act,1));
subplot(313);imagesc(mod(y-act-1,numStates)+1);
score = score(end);
%figure;hist(act(:),0:numStates*2);
%figure;plot(dpValues');hold all;plot(score);
