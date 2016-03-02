function [VValues,QValues,traj,acts, score] = spinnerOptimumP(fname,wavenums,VValues,QValues)
%use dynamic programming to determine each state's value and optimum action
%use probabilistic transitions - 1) where each action will lead on the next
%timestep (gaussian uncertainty in time-step); 2) which action an agent chooses (requires temperature param.)
%also differs from spinnerOptimum4 because it looks only n seconds into the
%future.

temperature = .1;%temperature for soft-max action selection (larger than .0015 explodes to inf)
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
dt = .1;
totallen = (wavelen + wavegap)*numel(wavenums);

travelcost = .0; %cost of moving per vertex, 0 appears to be OK!!
movecost = .01; %cost of moving at all -- prevents agent from taking a break while moving to a far target
popcost = 0; %cost of pop action -- discourages gratuitous pressing
spawnAccess = [poppabletime 1];
rampVal = linspace(spawnAccess(1),spawnAccess(2),diff(spawnAccess)*growtime/dt);

rawValues = zeros(numStates,totallen/dt);

for i = 1:size(dat,1)
    inds = round(dat(i,1)/dt + (1:numel(rampVal)) + poppabletime*growtime/dt);
    rawValues(dat(i,2)+1,inds) = rampVal;
end
if ~exist('VValues','var')
VValues = zeros(numStates,size(rawValues,2));
QValues = zeros(numStates,numStates*2,size(rawValues,2));
%isPopped = cell(size(VValues));
%act = zeros(size(rawValues));

uncertainT = .1; %std / sec;
thresh = .001;
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
for i = 1:size(rawValues,2)
    %take chunk of stimuli, blurred with temporal uncertainty
    rawChunk(:) = 0;
    visBubbles = find(round((dat(:,1)+growtime)/dt) >= i & round(dat(:,1)/dt) <= i);%%dat(:,1) > i*dt - growtime & dat(:,1) <= i*dt);
    if numel(visBubbles)
        if i >= 259
            3;
        end
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
                for k = 1:numStates
                    indsX(k) = mod(j-k-1,numStates)+1;
                    indsY([k k+numStates]) = min(size(Vlocal,2),round(ii+1+slidertime/dt*min(k,numStates-k)));
                    %indsY(k+numStates) = min(size(Vlocal,2),round(indsY(k)+popTime/dt)); %IF POPTIME IS NONZERO
                    %indsY([k k+numStates]) = min(size(Vlocal,2),indsY([k k+numStates]));
                    VBtemp(k,:) = VBlocal(indsX(k),:,indsY(k)); %Needed??
                    VBtemp(k+numStates,:) = VBtemp(k,:);
                    valNext([k k+numStates]) = Vlocal(indsX(k),indsY(k))-travelcost*min(k,numStates-k) - movecost*(k ~= numStates); %SINCE POPTIME = 0 NOW
                     if numel(curBubble) && stimChunk(curBubble,ii) > 0
                        VBtemp(k+numStates,curBubble) = stimChunk(curBubble,ii);
                        %valNext(k+numStates) = Vlocal(indsX(k),indsY(k+numStates))-travelcost*min(k,numStates-k)- movecost*(k ~= numStates) ...
                        valNext(k+numStates) = valNext(k+numStates) - popcost + stimChunk(curBubble,ii) - VBlocal(indsX(k),curBubble,indsY(k+numStates));
                    end
                end
                Qlocal(j,:,ii) = valNext;
                pAct = exp(valNext/temperature);
                if isinf(sum(pAct))
                    pAct = double(valNext == max(valNext));
                end
                pAct = pAct/sum(pAct);
                Vlocal(j,ii) = pAct*valNext';
                VBlocal(j,:,ii) = pAct*VBtemp;
            end
         end
%         if i >= 93
%             3
%         end
        VValues(:,i) = Vlocal(:,1);
        QValues(:,:,i) = Qlocal(:,:,1);
%        plot(Vlocal');title(i);drawnow;
    else
        QValues(:,1:numStates,i) = (eye(numStates)-1)*movecost;%ASSUMING NO TRAVELCOST
    end
end
end
traj = zeros(1,size(rawValues,2));
acts = traj;%acta = 1;
traj(1) = 1;
score = zeros(size(traj));
i = 2;%iOld = 1;
while i <= size(traj,2)
    %acta = act(traj(i-1),i-1);
%    if i == 
    [~,acta] = max(QValues(traj(i-1),:,i-1));%mod(acts(iOld)-1,numStates)+1
    acts(i) = acta;
    score(i) = score(i-1);
    if acta > numStates
        score(i) = score(i) + rawValues(traj(i-1),i-1);
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
%[~,y] = meshgrid(1:size(act,2),1:size(act,1));
%subplot(313);imagesc(mod(y-act-1,numStates)+1);
%figure;hist(act(:),0:numStates*2);
figure;plot(VValues');hold all;plot(score);
end

    
