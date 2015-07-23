function [Vd,Vm,runs] = embodiedValue(x_num, x_cost, t_cost)
%% plots an example of a value function and the associated bound, when you ...
%% have to move to the target
%
% sig2 is the overall task difficulty (variance of prior on mu), c is the
% evidence accumulation cost, and t is the time until which the values are
% to be computed. The function computes until 5*t that time, but only
% displays the value function / bounds until t.
%
% If not given, the arguments default to sig2 = 0.5^2, c = 0.1, and t = 3.
% 1) implement uncertainty about when the trial will be over... quick &
% dirty, have a supralinear cost associated with waiting too long.
% 2) implement uncertainty about where action places you in next time step
% 3) implement finer resolution x-position, to make derivatives (e.g.jerk)
% less bumpy
% 4) implement higher-dimensional state
% 5) implement costs that include higher-order costs (e.g. ^2, ^3)

%% settings
% task difficulty
if nargin < 1, x_num = 1; end%number of steps needed to move to either L or R target
% cost for moving
if nargin < 2, x_cost = 0.1; end;
% cost of time proceeding
if nargin < 3, t_cost = .1; end
% discretisation of belief and time (coarse, as only visualisation)
g_num = 50; %100
T = 1;
sig2 = .5^2;
dt = 0.0125;

%% compute the value function
makeFlat = @(x) x(:);

%% time steps
ts = 0:dt:T;
N = length(ts);

nInstances = 1000;
numAccum = 10;
x = meshgrid((1:N)-1,1:nInstances);
Vd = zeros(N,g_num,x_num*2+1);
ggOr = 0;Vhist = 0;
for i = 1:numAccum
    runs = randn(nInstances,length(ts)-1)*dt;
    rates = randn(nInstances,1)*sqrt(sig2)*dt;%(floor(rand(nInstances,1)*2)-.5)*dt/3;%
    runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
    runs = cumsum(runs,2);
    isRight = bsxfun(@times,ones(size(runs)),sign(rates))/2+.5;
    Vd(:,:,end) = squeeze(Vd(:,:,end)) + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],isRight(:),[N g_num],@sum);
    Vhist = Vhist + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],ones(1,numel(isRight)),[N g_num]);
    ggOr = ggOr + accumarray([makeFlat(x(:,1:end-1))+1 binFun(makeFlat(runs(:,1:end-1)),g_num) ...
        binFun(makeFlat(runs(:,2:end)),g_num)],ones(1,(N-1)*nInstances),[N-1 g_num g_num]);%...
        %min(g_num,max(1,round((makeFlat(runs(:,1:end-1))/scale+.5)*g_num))) ...
        %min(g_num,max(1,round((makeFlat(runs(:,2:end))/scale+.5)*g_num)))],ones(1,(N-1)*nInstances));
end
Vd(:,:,end) = squeeze(Vd(:,:,end))./Vhist;
%runs = binFun(runs,g_num);
%for i = 1:size(ggOr,2)
%     for j = 1:size(ggOr,3)
%         ggOr(:,i,j) = filtfilt(gausswin(7),sum(gausswin(7)),ggOr(:,i,j));
%
%end

Vd(isnan(Vd)) = 0.5;
for i = 1:size(Vd,2)
    Vd(:,i,end) = filtfilt(gausswin(.25/dt),sum(gausswin(.25/dt)),Vd(:,i,end));
end
Vd(:,:,1) = 1-Vd(:,:,end);
ggOr = bsxfun(@rdivide, ggOr, sum(ggOr, 3)+eps);
Vm = NaN(N-1, g_num,x_num*2+1);

for i = N-1:-1:1
    evidence = squeeze(ggOr(i,:,:))*squeeze(Vd(i+1,:,:));
    timeCosts = -t_cost*dt*[0 ones(1,2*x_num-1) 0];    
    for j = 2:2*x_num
        moveCosts = -x_cost*((1:2*x_num+1)-j).^2/(x_num).^2;%
        [Vd(i,:,j) Vm(i,:,j)] = max(evidence+repmat(moveCosts+timeCosts,[g_num 1]),[],2);
    end
end

function dat = binFun(dat,nBins)
dat = round((tanh(dat)+1)/2*(nBins-1)+1);