function [Vd,Vm,RTs,Rchoice] = embodiedValueMinimal1(x_num, y_num, sig2, x_cost, t_cost)
%% calculates value function for embodied agent to make a choice
% INPUTS
% x_num - number of steps it takes for agent to get to the target [1]
% y_num - maximum number of time steps that trial can last [10]
% sig2 - overall task difficulty (variance of prior on mu) [.5^2]
% x_cost - cost of pressing a button [0]
% t_cost - evidence accumulation cost per dt [1/y_num]

%% settings
%number of steps needed to move to either L or R target
if nargin < 1, x_num = 1; end
%number of steps needed to get to end of trial
if nargin < 2, y_num = 10; end 
% variance of prior on mu
if nargin < 3, sig2 = 1^2; end
% cost for moving (i.e. button press)
if nargin < 4, x_cost = 0; end;
% cost of time proceeding
if nargin < 5, t_cost = 1/y_num; end
% discretisation of belief and time (coarse, as only visualisation)
g_num = 100;

dt = .0125;
T = dt*y_num;
ts = 0:dt:T;
N = length(ts);

dg = linspace(1 / g_num / 2, 1 - 1 / g_num / 2, g_num);
invgs = norminv(dg);

Vd = zeros(N,g_num,x_num*2+1);
Vd(:,:,end) = repmat(dg,[N 1]);
Vd(:,:,1) = 1-Vd(:,:,end);
%Vd(:,:,2) = max(Vd(:,:,1),Vd(:,:,end)); %%REMOVE THIS
Vm = NaN(N-1, g_num,x_num*2+1);

for i = N-1:-1:1
    gg = belieftrans(invgs, dt / (ts(i) + 1/sig2));
    evidence = gg*squeeze(Vd(i+1,:,:));
    timeCosts = -t_cost*[0 ones(1,2*x_num-1) 0];  
    for j = 2:2*x_num
        moveCosts = abs((1:2*x_num+1)-j);
        moveCosts(moveCosts > 1) = Inf;
        moveCosts = -x_cost*moveCosts;
        [Vd(i,:,j), Vm(i,:,j)] = max(evidence+repmat(moveCosts+timeCosts,[g_num 1])+randn(size(evidence))*eps,[],2);
    end
end

%% make some example runs
nInstances = 200;
runs = randn(nInstances,length(ts)-1)*sqrt(dt);
rates = randn(nInstances,1)*sqrt(sig2)*dt;
runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
runs = cumsum(runs,2);
p = normcdf(bsxfun(@rdivide,runs,sqrt(ts + 1/sig2)));
runsBin = round((g_num-1)*p)+1;
pos = nan(size(runs,1),size(runs,2));
pos(:,1) = x_num+1;
for i = 2:size(runsBin,2)
    for j = 1:size(runsBin,1)
        if ~isnan(pos(j,i-1))
            pos(j,i) = Vm(i-1,runsBin(j,i-1),pos(j,i-1));
        end
    end
end
%bounds = [-inf,2,-1,0,1,2,inf];%-2,-1,0,1,2
c = jet(6);
prc = linspace(0,100,size(c,1)+1);
cc1 = zeros(nInstances,3);
group = zeros(nInstances,1);
for i = 1:size(c,1)
inds = rates >= prctile(rates,prc(i)) & rates <= prctile(rates,prc(i+1));
group(inds) = i;
%rates/dt/sqrt(sig2) > bounds(i) & rates/dt/sqrt(sig2) < bounds(i+1);
cc1(inds,:) = repmat(c(i,:),[sum(inds) 1]);
%subplot(2,3,i);
%plot((pos(inds,:)+randn(size(pos(inds,:)))*.1)',ts,'color',c(i,:));hold all;
%plot(nanmean(pos(inds,:)),ts,'color',c(i,:));hold all;
%corr(i) = mean(sign(runs(inds,end)) > 0); 
%nums(i) = sum(inds);
end
runs = zeros(nInstances,2);
r = 1:nInstances;%randperm(nInstances);
%lastPos = zeros(nInstances,1);
for i= 1:nInstances
    ind = find(~isnan(pos(i,:)),1,'last');
    runs(i,1) = pos(i,ind) > x_num;
    runs(i,2) = ind*dt;
%    plot((pos(r(i),:)-x_num-1)/x_num+randn(1)*.01,ts+randn(1)*.01*dt,'color',cc1(r(i),:),'linewidth',1);hold all;
end
%axis tight;
%set(gca,'xlim',[-1 1]);%,'ylim',[0 max(ts)]);
RTs = accumarray(group,runs(:,2),[],@mean);
Rchoice = accumarray(group,runs(:,1),[],@mean);
%figure;subplot(211);plot(RTs);subplot(212);plot(Rchoice);


function gg = belieftrans(invgs, dteff)
%% Returns the belief transition matrix p(g' | g, t)
invgdiff = bsxfun(@minus, invgs, sqrt(1 + dteff) * invgs');
gg = exp(bsxfun(@minus, invgs.^2 / 2, invgdiff.^2 / (2 * dteff)));
gg = bsxfun(@rdivide, gg, sum(gg, 2));