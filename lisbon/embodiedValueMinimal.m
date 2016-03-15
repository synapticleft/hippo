function [Vd,Vm,runsBin,p] = embodiedValueMinimal(x_num, y_num, x_cost, t_cost)
%% plots an example of a value function and the associated bound, when you ...
%% have to move to the target
%
% sig2 is the overall task difficulty (variance of prior on mu), t_cost is the
% evidence accumulation cost, and t is the time until which the values are
% to be computed.

% If not given, the arguments default to sig2 = 0.5^2, x_num = 1, y_num = 0,
% t_cost = 0.1.
% how should choice value depend on time elapsed?
% - (1 or 0) - timeElapsed
% - (1 - timeElapsed or -(1 - timeElapsed)
% - (1 - timeElapsed or -1)

%% settings
if nargin < 1, x_num = 1; end %number of steps needed to move to either L or R target
if nargin < 2, y_num = 10; end %number of steps needed to get to end of trial
% cost for moving
if nargin < 3, x_cost = eps; end;
% cost of time proceeding
if nargin < 4, t_cost = 1/y_num; end
% discretisation of belief and time (coarse, as only visualisation)
g_num = 51; %100
sig2 = 1^2;
dt = .0125;
T = dt*y_num;

%% time steps
ts = 0:dt:T;
dg = linspace(1 / g_num / 2, 1 - 1 / g_num / 2, g_num);
invgs = norminv(dg);
%dg = 1/(g_num*2):1/g_num:1-1/(g_num*2);%0:1/g_num:1;%0:1/g_num:1;%
N = length(ts);
Vd = zeros(N,g_num,x_num*2+1);
   
Vd(:,:,end) = repmat(dg,[N 1]);
Vd(:,:,1) = 1-Vd(:,:,end);
%ggOr = bsxfun(@rdivide, ggOr, sum(ggOr, 3)+eps);
Vm = NaN(N-1, g_num,x_num*2+1);

for i = N-1:-1:1
    gg = belieftrans(invgs, dt / (ts(i) + sig2));
    evidence = gg*squeeze(Vd(i+1,:,:));%squeeze(ggOr(i,:,:))
    timeCosts = -t_cost*[0 ones(1,2*x_num-1) 0];%dt*    
    for j = 2:2*x_num
        moveCosts = abs((1:2*x_num+1)-j);
        moveCosts(moveCosts > 1) = Inf;
        moveCosts = -x_cost*moveCosts;%((1:2*x_num+1)-j).^2/(x_num).^2;%
        [Vd(i,:,j), Vm(i,:,j)] = max(evidence+repmat(moveCosts+timeCosts,[g_num 1]),[],2);
    end
end
%% make some example runs
nInstances = 100;
runs = randn(nInstances,length(ts)-1)*sqrt(dt);
rates = randn(nInstances,1)*sqrt(sig2)*dt;%(floor(rand(nInstances,1)*2)-.5)*dt/3;%
runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
runs = cumsum(runs,2);
%p = 1-normcdf(0, runs./(1/sig2+dt*repmat(ts,size(runs,1),1)),1./(1/sig2+dt*repmat(ts,size(runs,1),1)));
p = normcdf(bsxfun(@rdivide,runs,sqrt(ts + 1/sig2)));
runsBin = round((g_num-1)*p)+1;
%runsBin = runs;%binFun(runs,g_num);%round((tanh(runs)+1)/2*(50-1)+1);
pos = nan(size(runs,1),size(runs,2));
pos(:,1) = x_num+1;
for i = 2:size(runsBin,2)
    for j = 1:size(runsBin,1)
        if ~isnan(pos(j,i-1))
            pos(j,i) = Vm(i-1,runsBin(j,i-1),pos(j,i-1));
        end
    end
end
%figure;subplot(211);plot(runsBin');
%subplot(212);plot(pos');
plot(pos',ts);set(gca,'xlim',[1 x_num*2+1],'ylim',[0 max(ts)]);

function gg = belieftrans(invgs, dteff)
%% Returns the belief transition matrix p(g' | g, t)
invgdiff = bsxfun(@minus, invgs, sqrt(1 + dteff) * invgs');
gg = exp(bsxfun(@minus, invgs.^2 / 2, invgdiff.^2 / (2 * dteff)));
gg = bsxfun(@rdivide, gg, sum(gg, 2));


    
% %%OLD
% ggOr = zeros(N,g_num,g_num);
% ggOr1 = ggOr;
% for i = 1:numel(ts)
%     for g = 1:g_num 
%         ggOr(i,g,:)= BeliefTransitionDrugo(dg,ts(i),dt,dg(g),1/sig2);
%     end
% %    ggOr(i,:,:) = belieftrans(invgs, dt / (ts(i) + sig2));
% end 
% function gg = BeliefTransitionDrugo( gj, t, tao, gi,sigma2 )
% invgs = norminv(gi);
% invgsp = norminv(gj);
% Steff = tao ./ (t + 1/sigma2);
% invgdiff = invgsp - sqrt(1 + Steff) .* invgs;
% % unnormalised gg
% ggu = 1./sqrt(Steff)*exp( invgsp.^2 / 2 - invgdiff.^2 ./ (2 * Steff));
% gg = ggu./sum(ggu,2);

%function dat = binFun(dat,nBins)
%dat = round((tanh(dat)+1)/2*(nBins-1)+1);

% gausSmooth = 2;
%% compute the value function
% makeFlat = @(x) x(:);
% nInstances = 1000;
% numAccum = 1000;
% x = meshgrid((1:N)-1,1:nInstances);
% ggOr = 0;Vhist = 0;
% for i = 1:numAccum
%     runs = randn(nInstances,length(ts)-1)*sqrt(dt);
%     rates = randn(nInstances,1)*sqrt(sig2)*dt;%(floor(rand(nInstances,1)*2)-.5)*dt/3;%
%     runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
%     runs = binFun(cumsum(runs,2),g_num);
%     isRight = bsxfun(@times,ones(size(runs)),sign(rates))/2+.5;
%     Vd(:,:,end) = squeeze(Vd(:,:,end)) + accumarray([x(:)+1 runs(:)],isRight(:),[N g_num],@sum);
%     Vhist = Vhist + accumarray([x(:)+1 runs(:)],ones(1,numel(isRight)),[N g_num]);
%     ggOr = ggOr + accumarray([makeFlat(x(:,1:end-1))+1 makeFlat(runs(:,1:end-1)) ...
%         makeFlat(runs(:,2:end))],ones(1,(N-1)*nInstances),[N-1 g_num g_num]);%...
%         %min(g_num,max(1,round((makeFlat(runs(:,1:end-1))/scale+.5)*g_num))) ...
%         %min(g_num,max(1,round((makeFlat(runs(:,2:end))/scale+.5)*g_num)))],ones(1,(N-1)*nInstances));
% end
% % for i = 1:size(Vd,2)
% %     Vd(:,i,end) = filtfilt(gausswin(gausSmooth),sum(gausswin(gausSmooth)),Vd(:,i,end));
% %     Vhist(:,i) = filtfilt(gausswin(gausSmooth),sum(gausswin(gausSmooth)),Vhist(:,i));
% % end
% % for i = 1:size(Vd,1)
% %     Vd(i,:,end) = filtfilt(gausswin(gausSmooth),sum(gausswin(gausSmooth)),Vd(i,:,end));
% %     Vhist(i,:) = filtfilt(gausswin(gausSmooth),sum(gausswin(gausSmooth)),Vhist(i,:));
% % end
% Vd(:,:,end) = squeeze(Vd(:,:,end))./Vhist;
% %runs = binFun(runs,g_num);
% %for i = 1:size(ggOr,2)
% %     for j = 1:size(ggOr,3)
% %         ggOr(:,i,j) = filtfilt(gausswin(7),sum(gausswin(75)),ggOr(:,i,j));
% %
% %end
% 
% Vd(isnan(Vd)) = 0.5;