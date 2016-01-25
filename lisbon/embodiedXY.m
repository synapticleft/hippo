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
% less bumpy -- done!
% 4) implement higher-dimensional state -- done!
% 5) implement costs that include higher-order costs (e.g. ^2, ^3)
% 6) implement y- dimension

makeFlat = @(x) x(:);
binFun = @(dat,nBins) round((tanh(dat)+1)/2*(nBins-1)+1);

    
%% settings
% task difficulty
x_num = 3; %number of steps needed to move to either L or R target
y_num = 3; 

% cost for moving
dist_cost = 1; 

% cost of time proceeding
t_cost = .1; 
% discretisation of belief and time (coarse, as only visualisation)
g_num = 50; %100
T = 3;
sig2 = 1^2;
dt = 0.0125;
ord = 1;

%% compute the value function
makeFlat = @(x) x(:);

%% time steps
ts = 0:dt:T;
N = length(ts);

nInstances = 1000;
numAccum = 500;
x = meshgrid((1:N)-1,1:nInstances);
Vd = zeros(N,g_num,x_num*2+1,y_num);

ggOr = 0;Vhist = 0;
for i = 1:numAccum
%     runs = randn(nInstances,length(ts)-1)*sqrt(dt);
%     rates = randn(nInstances,1)*sqrt(sig2)*dt;%(floor(rand(nInstances,1)*2)-.5)*dt/3;%
%     runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
    runs = randn(nInstances,length(ts)-1)*sqrt(sig2*dt);
    rates = randn(nInstances,1)*dt;
    runs = [zeros(nInstances,1) bsxfun(@plus,runs,rates)];
    runs = cumsum(runs,2);
    isRight = bsxfun(@times,ones(size(runs)),sign(rates))/2+.5;
    %R: V(t,g,xpos,ypos)
    Vd(:,:,end,end) = squeeze(Vd(:,:,end,end)) + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],isRight(:),[N g_num],@sum);
    Vhist = Vhist + accumarray([x(:)+1 binFun(makeFlat(runs(:)),g_num)],ones(1,numel(isRight)),[N g_num]);
    ggOr = ggOr + accumarray([makeFlat(x(:,1:end-1))+1 binFun(makeFlat(runs(:,1:end-1)),g_num) ...
        binFun(makeFlat(runs(:,2:end)),g_num)],ones(1,(N-1)*nInstances),[N-1 g_num g_num]);%...
        %min(g_num,max(1,round((makeFlat(runs(:,1:end-1))/scale+.5)*g_num))) ...
        %min(g_num,max(1,round((makeFlat(runs(:,2:end))/scale+.5)*g_num)))],ones(1,(N-1)*nInstances));
end
Vd(:,:,end,end) = squeeze(Vd(:,:,end,end))./Vhist;
%runs = binFun(runs,g_num);
%for i = 1:size(ggOr,2)
%     for j = 1:size(ggOr,3)
%         ggOr(:,i,j) = filtfilt(gausswin(7),sum(gausswin(7)),ggOr(:,i,j));
%
%end

Vd(isnan(Vd)) = 0.5;
for i = 1:size(Vd,2)
    Vd(:,i,end,end) = filtfilt(gausswin(.25/dt),sum(gausswin(.25/dt)),Vd(:,i,end,end));
end
%% should only still endings be rewarded, or is it OK to keep moving?
%% I think only still, since other values would require an overshoot in x (which is bounded)
%Vd(:,:,1:end-1,end) = repmat(Vd(:,:,end,end),[1 1 size(Vd,3)-1]);
%R: Mirror values for xpos 1 (left), only for terminal y 
Vd(:,:,1,end) = 1-Vd(:,:,end,end);
ggOr = bsxfun(@rdivide, ggOr, sum(ggOr, 3)+eps);
Vm = NaN(N-1, g_num,x_num*2+1,y_num);
%Vdd = Vm;
% Vdd = repmat(Vd,[1 1 1 1 size(Vd,3) size(Vd,4)]);
Vdd = repmat(Vd,[1 1 1 1 ]);

for i = N-1:-1:1 %back propagation in time (i)
    gridVal = ones(2*x_num+1,y_num);
    gridVal([1, 2*x_num+1],y_num)=[0, 0];
    timeCosts = -t_cost*dt*gridVal;  %R: time costs for each xy in grid
    
        for k = 1:1:y_num%R: At each time, compute evidence (belief) for all x and all y (k)
            %evidence = squeeze(ggOr(i,:,:))*squeeze(Vdd(i+1,:,:,k));
            evidence(:,:,k) = squeeze(ggOr(i,:,:))*squeeze(Vdd(i+1,:,:,k));
        end
        
        for k = y_num:-1:1  %R: back propagation in time running over all x (j) and y (k)
            for j = 1:2*x_num+1
                
            for y=1:y_num %cost for jump will be proportional to euclidean distance traveled (normalized by grid size)
                 for x=1:2*x_num+1
                    distance(x,y) = sqrt((((y-k)/y_num)^2)+(((x-j)/(2*x_num+1))^2));
                 end      
            end
                moveCosts = -dist_cost(1)*distance.^2;
                
            for g = 1:g_num    
                ExpectedValues = squeeze(evidence(g,:,:))+moveCosts+timeCosts;  
                if (j == 1 && k == y_num) || (j==2*x_num+1 && k == y_num)
                    %do not overwrite value of end states               
                else
                [Vdd(i,g,j,k) Vm(i,g,j,k)] = max(ExpectedValues(:));
                end
            end
            end
        end
end
% 

gs = [1 10 26 40 50];
figure
for t = size(Vdd,1):-1:1
    ind = 1;
    for g = gs
    subplot(1,length(gs),ind)
    vals = squeeze(Vdd(t,g,:,:));
    imagesc(vals')
    xlabel('x')
    ylabel('y')
    if g==g_num
        title(strcat('g=', num2str(1-1/g_num)))
    else
        title(strcat('g=',num2str(g/(g_num))))
    end
    ind = ind+1;
    end
    if t==size(Vdd,1)
      pause(2)
    else 
      pause(.05)
    end
    
end

% gs = [1 10 26 40 50];
% 
% movieObject=VideoWriter(strcat('SimpleAgentXY','.avi'));
% movieObject.FrameRate=30;
% open(movieObject)
% 
% figure
% for t = size(Vdd,1):-1:1
%     ind = 1;
%     for g = gs
%     subplot(1,length(gs),ind)
%     vals = squeeze(Vdd(t,g,:,:));
%     imagesc(vals')
%     xlabel('x')
%     ylabel('y')
%     if g==g_num
%         title(strcat('g=', num2str(1-1/g_num)))
%     else
%         title(strcat('g=',num2str(g/(g_num))))
%     end
%     ind = ind+1;
%     end
%     if t==size(Vdd,1)
%       pause(2)
%     else 
%       pause(.05)
%     end
%     
%         
%     movieframe=getframe(gcf);
%     writeVideo(movieObject, movieframe)
% end
% close(movieObject);
% 

% 
% figure
% for t = 1:size(Vdd,1)
%     for y = 1:size(Vdd,4)
%     subplot(1,y_num,y)
%     vals = squeeze(Vdd(t,:,:,y));
%     imagesc(vals)
%     end
% waitforbuttonpress
% end
% 
% 
% 
% figure
% for t = 1:size(Vdd,1)
%     for x = 1:size(Vdd,3)
%     subplot(1,2*x_num+1,x)
%     vals = squeeze(Vdd(t,:,x,:));
%     imagesc(vals)
%     end
% waitforbuttonpress
% end
% %% make some example runs
% if exist('runs1','var')
%     runsBin = runs1;
% else
%     runsBin = binFun(runs,g_num);%round((tanh(runs)+1)/2*(50-1)+1);
% end
% p = nan*ones(size(runs,1),size(runs,2),ord+1);
% p(:,1,:) = x_num+1;
% for i = 2:size(runsBin,2)
%     for j = 1:size(runsBin,1)
%         if ~isnan(p(j,i-1,end)) %&& ~isnan(p(j,i-1,1))
%             p(j,i,1:end-1) = p(j,i-1,2:end);
%             p(j,i,end) = Vm(i-1,runsBin(j,i-1),p(j,i-1,2),p(j,i-1,1));
%         end
%     end
% end
% 