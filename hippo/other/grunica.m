function [w wz] = grunica(x)

Ls = .0001;%[.01 .001 .0001];
sweeps = 1000;
x = bsxfun(@rdivide,x,std(x,0,2));
[N P] = size(x);
wz= 2*pinv(sqrtm(cov(x')));                  % get decorrelating matrix
x=wz*x;                              % decorrelate mixes so cov(x')=4*eye(N);
Id=eye(N);
w= eye(N);%rand(N);                            % init. unmixing matrix, or w=rand(M,N);
oldw=w; olddelta=ones(1,N*N);
lrate=0.00065/log(N); %0.0001; 
B=min(floor(5*log(P)),floor(0.3*P))%2*64; 
sweeps = 300;
annealdeg = 60;annealstep = 0.90;degconst = 180/pi;nochange = 0.0000001; 
change = eps;oldchange  = change;

for ii = 1:sweeps %numel(Ls)
x=x(:,randperm(P));
noblocks=fix(P/B);
BI=B*Id;
for t=1:B:noblocks*B,
  u=w*x(:,t:t+B-1); 
  w=w+lrate*(BI+(1-2*(1./(1+exp(-u))))*u')*w;
end
      delta=w(:).'-oldw(:).';
      change=delta*delta'; 
  %[change,oldw,olddelta,angle] = sepout(w,oldw,olddelta,lrate);
  angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
  fprintf('step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n',sweeps,lrate,change,degconst*angledelta);
  oldw = w;
      if degconst*angledelta > annealdeg,  
        lrate = lrate*annealstep;          % anneal learning rate
        olddelta   = delta;                % accumulate angledelta until
        oldchange  = change;               %  annealdeg is reached
      end
       if sweeps >2 && change < nochange,      % apply stopping rule
           break
%         laststep=step;            
%            step=maxsteps;                  % stop when weights stabilize
%       elseif change > DEFAULT_BLOWUP,      % if weights blow up,
%         lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
       end                                 % with a smaller learning rate
end
%w = w*wz;

function [change,oldw,olddelta,angle] = sepout(w,oldw,olddelta,lrate)
[change,olddelta,angle]=wchange(oldw,w,olddelta); 
oldw=w;
fprintf('****lrate=%.6f change=%.6f angle=%.2f deg. \n',...
   lrate,change,180*angle/pi);

function [change,delta,angle]=wchange(w,oldw,olddelta)
  [M,N]=size(w); delta=reshape(oldw-w,1,M*N);
  change=delta*delta';
  angle=acos((delta*olddelta')/sqrt((delta*delta')*(olddelta*olddelta')));