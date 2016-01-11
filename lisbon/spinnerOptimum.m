function spinnerOptimum(fname,wavenum)
%use dynamic programming to determine each state's value and optimum action


dat = importdata(fname);
dat = dat(dat(:,3) == wavenum,1:2);
numStates = numel(unique(dat(:,2)));
dt = .01;
growtime = 5;
wavelen = 30 + growtime;
slidertime = .4;
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
    %inds = spawnAccess(1)*growConst:spawnAccess(2)*growConst;
%    for j = 1:numStates
%        inds = round(abs(dat(i,2)+1-j)*slidertime/dt + dat(i,1)/dt + (1:numel(rampVal)));
%        inds(inds > size(rawValues,2))= [];
%        rawValues(j,inds) = rawValues(j,inds) + rampVal(1:numel(inds)); 
%%        rawValues(dat(i,2)+1,dat(i,1)/dt+inds) = ;
%    end
%    rawValues(dat(i,2)+1,dat(i,1)/dt+inds) = linspace(spawnAccess(1),spawnAccess(2),numel(inds));
end
dpValues = zeros(size(rawValues,1),size(rawValues,2)+1);
isPopped = cell(size(rawValues,1),size(rawValues,2)+1);
%act(:,end) = 1:numStates;dpValues = rawValues(:,end);
%act(:,end) = 0;dpValues(:,end) = 0;
numActs= 4;
for i = size(rawValues,2):-1:1
    for j = 1:numStates
        valNext(1) = dpValues(mod(j-1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))) - travelcost;
        valNext(2) = dpValues(mod(j+1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))) - travelcost;
        valNext(3) = dpValues(j,i+1);
        if ~ismember(theBubble(j,i),isPopped{j,i+1})
            valNext(4) = rawValues(j,i) + valNext(3);
            flag = 3;
        elseif ~ismember(theBubble(j,i),isPopped{mod(j-1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))})
            valNext(4) = rawValues(j,i) + valNext(1);
            flag = 1;
        elseif ~ismember(theBubble(j,i),isPopped{mod(j+1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))})
            valNext(4) = rawValues(j,i) + valNext(2);
            flag = 2;
        else
            valNext(4) = -inf;
        end
        %for k = 1:numStates
        %valNext(k) = dpValues(k,round(min(size(dpValues,2),i+1+abs(k-j)*slidertime/dt))) - abs(k-j)*travelcost;
        %end
        %valNext = valNext + rawValues(j,i); %
        %valNext(j) = valNext(j) - rawValues(j,i);
        [dpValues(j,i),act(j,i)] = max(valNext);
        if act(j,i) == 4
            if flag == 3
                isPopped{j,i} = union(theBubble(j,i),isPopped{j,i+1});
            elseif flag == 1
                isPopped{j,i} = union(theBubble(j,i),isPopped{mod(j-1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))});
            elseif flag == 2
                isPopped{j,i} = union(theBubble(j,i),isPopped{mod(j+1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))});
            end
        elseif act(j,i) == 3
            isPopped{j,i} = isPopped{j,i+1};
        elseif act(j,i) == 2
            isPopped{j,i} = isPopped{mod(j+1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))};
        elseif act(j,i) == 1
            isPopped{j,i} = isPopped{mod(j-1-1,numStates)+1,round(min(size(dpValues,2),i+1+slidertime/dt))};
        end
    end
end
figure;subplot(311);imagesc(rawValues);
subplot(312);imagesc(dpValues);
subplot(313);imagesc(act);
end

