function [Vd,Ve] = embodiedValue(sig2, c, T)
%% plots an example of a value function and the associated bound, when you ...
%% have to move to the target
%
% sig2 is the overall task difficulty (variance of prior on mu), c is the
% evidence accumulation cost, and t is the time until which the values are
% to be computed. The function computes until 5*t that time, but only
% displays the value function / bounds until t.
%
% If not given, the arguments default to sig2 = 0.5^2, c = 0.1, and t = 3.


%% settings
% task difficulty
if nargin < 1, sig2 = 0.5^2; end;
% cost for accumulating evidence
if nargin < 2, c = 0.1; end;
% time-frame of interest
if nargin < 3, T = 3; end
% discretisation of belief and time (coarse, as only visualisation)
g_num = 50; %100
x_num = 5; %number of steps needed to move to either L or R target
dt = 0.1;

%% compute the value function
makeFlat = @(x) x(:);

%% time steps
ts = 0:dt:T;
N = length(ts);

% ggAll = zeros(numel(gs),numel(gs),N);
% %bins = f;
nInstances = 10000;
numAccum = 10;
x = meshgrid((1:N)-1,1:nInstances);
Vd = 0;
ggOr = 0;Vhist = 0;
for i = 1:numAccum
    runs = randn(nInstances,length(ts)-1)*dt;
    rates = randn(nInstances,1)*sqrt(sig2)*dt;%(floor(rand(nInstances,1)*2)-.5)*dt/3;%
    runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
    runs = cumsum(runs,2);
    isRight = bsxfun(@times,sign(runs),sign(rates))/2+.5;
    Vd = Vd + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],isRight(:),[],@sum);
    Vhist = Vhist + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],ones(1,numel(isRight)));
    ggOr = ggOr + accumarray([makeFlat(x(:,1:end-1))+1 binFun(makeFlat(runs(:,1:end-1)),g_num) ...
        binFun(makeFlat(runs(:,2:end)),g_num)],ones(1,(N-1)*nInstances));%...
        %min(g_num,max(1,round((makeFlat(runs(:,1:end-1))/scale+.5)*g_num))) ...
        %min(g_num,max(1,round((makeFlat(runs(:,2:end))/scale+.5)*g_num)))],ones(1,(N-1)*nInstances));
end
Vd = Vd./Vhist;Vd(isnan(Vd)) = 1;

%for i = 1:size(ggOr,2)
%     for j = 1:size(ggOr,3)
%         ggOr(:,i,j) = filtfilt(gausswin(7),sum(gausswin(7)),ggOr(:,i,j));
%
%end
for i = 1:size(Vd,2)
    Vd(:,i) = filtfilt(gausswin(1/dt),sum(gausswin(1/dt)),Vd(:,i));
end
ggOr = bsxfun(@rdivide, ggOr, sum(ggOr, 3));
figure;imagesc(Vhist);
 figure;imagesc(Vd);
% figure;imagesc(squeeze(ggOr(end,:,:))-squeeze(ggOr(end,:,:))');
% figure;imagesc(squeeze(ggOr(round(end*.75),:,:) - ggOr(round(end/4))));
% figure;imagesc(ggAll-ggAll');
% return
Ve = NaN(N-1, g_num);
Vd = Vd(2:end,:,:);
for i = N-1:-1:1
    % compute next value, V(g, t(i) + dt)
    if i == N-1
        % in the last step, there is no next Ve
        Vt = Vd(i,:);
%        Vt1 = Vd(i,:);
    else
        Vt = max(Vd(i,:), Ve(i+1, :));
%        Vt1 = max(Vd(i,:),Ve1(i+1,:));
    end
    % based on this, compute current value for accumulating more evidence
    gg = squeeze(ggOr(i,:,:));%belieftrans(invgs, dt / (ts(i) + sig2));
    % Vt1 * gg' results in the expected future value
    Ve(i, :) = Vt * gg' - c * dt;
%    Ve1(i,:) = Vt1* ggAll' - c * dt;
end

function dat = binFun(dat,nBins)
dat = round((tanh(dat)+1)/2*(nBins-1)+1);