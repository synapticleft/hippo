function [Vd,Vm,runsBin,p] = embodiedValueMinimal1(x_num, y_num, sig2, x_cost, t_cost)
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
if nargin < 3, sig2 = .5^2; end
% cost for moving (i.e. button press)
if nargin < 4, x_cost = 0; end;
% cost of time proceeding
if nargin < 5, t_cost = 1/y_num; end
% discretisation of belief and time (coarse, as only visualisation)
g_num = 51;

dt = .0125;
T = dt*y_num;
ts = 0:dt:T;
N = length(ts);

dg = linspace(1 / g_num / 2, 1 - 1 / g_num / 2, g_num);
invgs = norminv(dg);

Vd = zeros(N,g_num,x_num*2+1);
Vd(:,:,end) = repmat(dg,[N 1]);
Vd(:,:,1) = 1-Vd(:,:,end);
Vm = NaN(N-1, g_num,x_num*2+1);

for i = N-1:-1:1
    gg = belieftrans(invgs, dt / (ts(i) + sig2));
    evidence = gg*squeeze(Vd(i+1,:,:));
    timeCosts = -t_cost*[0 ones(1,2*x_num-1) 0];  
    for j = 2:2*x_num
        moveCosts = abs((1:2*x_num+1)-j);
        moveCosts(moveCosts > 1) = Inf;
        moveCosts = -x_cost*moveCosts;
        [Vd(i,:,j), Vm(i,:,j)] = max(evidence+repmat(moveCosts+timeCosts,[g_num 1]),[],2);
    end
end

%% make some example runs
nInstances = 100;
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
plot(pos',ts);set(gca,'xlim',[1 x_num*2+1],'ylim',[0 max(ts)]);

function gg = belieftrans(invgs, dteff)
%% Returns the belief transition matrix p(g' | g, t)
invgdiff = bsxfun(@minus, invgs, sqrt(1 + dteff) * invgs');
gg = exp(bsxfun(@minus, invgs.^2 / 2, invgdiff.^2 / (2 * dteff)));
gg = bsxfun(@rdivide, gg, sum(gg, 2));