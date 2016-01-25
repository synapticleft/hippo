function spinnerOptimum1(fname,wavenums)
%use dynamic programming to determine each state's value and optimum action

%r =  [.1 .2 .3 .4 .5 .4 .3 .2];%[.1 .3 .1 .2 .1 .4 .5 .6];

dat = importdata(fname);
dat = dat(ismember(dat(:,3), wavenums),1:2);
%dat = [];dat(:,1) = ones(1,8)+r;dat(:,2) = 0:7;%+rand(1,8)
numStates = numel(unique(dat(:,2)));
dt = .01;
growtime = 5;
wavelen = 30 + growtime;
slidertime = .45;
travelcost = .01;
spawnAccess = [.35 1];

rawValues = zeros(numStates,wavelen/dt);
growConst = growtime/dt;
rampVal = linspace(spawnAccess(1),spawnAccess(2),diff(spawnAccess)*growConst);
theBubble = rawValues;
for i = 1:size(dat,1)
    inds = round(dat(i,1)/dt + (1:numel(rampVal)));
    rawValues(dat(i,2)+1,inds) = rampVal;
    theBubble(dat(i,2)+1,inds) = i;
end
dpValues = zeros(size(rawValues,1),size(rawValues,2)+1);
isPopped = cell(size(rawValues,1),size(rawValues,2)+1);
act = zeros(size(rawValues));
for i = size(rawValues,2):-1:1
    for j = 1:numStates
        indsX = [mod(j-1-1,numStates)+1 mod(j+1-1,numStates)+1 j];
        indsY = [round(min(size(dpValues,2),i+1+slidertime/dt)) round(min(size(dpValues,2),i+1+slidertime/dt)) i+1];
        for k = 1:3
            valNext(k) = dpValues(indsX(k),indsY(k));
            if k < 2
                valNext(k) = valNext(k) - travelcost;
            end
            if ~ismember(theBubble(j,i),isPopped{indsX(k),indsY(k)})
                valNext(k+3) = valNext(k) + rawValues(j,i);
            else
                valNext(k+3) = -Inf;
            end
        end
        [dpValues(j,i),act(j,i)] = max(valNext);
        if act(j,i) >= 4
            isPopped{j,i} = union(theBubble(j,i),isPopped{indsX(act(j,i)-3),indsY(act(j,i)-3)});
        else
            isPopped{j,i} = isPopped{indsX(act(j,i)),indsY(act(j,i))};
        end
    end
end
traj = zeros(1,size(rawValues,2));
traj(1) = 1;
score = zeros(size(traj));
sh = [-1 1 0];
for i = 2:size(traj,2)
    traj(i) = mod(traj(i-1)+sh(mod(act(traj(i-1),i-1)-1,3)+1)-1,numStates)+1;
    score(i) = score(i-1);
    if act(traj(i-1),i-1) >= 4
        score(i) = score(i) + rawValues(traj(i),i);
    end
end
figure;plot(score);
figure;plot(traj);
% %for i = 1:size(traj,1)
% %    for j = 1:size(traj,2)
% %        if j == 1
% %            traj(i,j) = 1;
% %        end
% %    end
% %end
% figure;subplot(311);imagesc(rawValues);
% subplot(312);imagesc(dpValues);
% subplot(313);imagesc(act);
% figure;hist(act(:),0:6);
% figure;plot(dpValues');
end

