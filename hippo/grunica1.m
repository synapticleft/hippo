function [weights,sphere,activations] = grunica1(data,sphering)
%% complex valued infomax, distilled

[chans frames] = size(data); % determine the data size
DEFAULT_BLOWUP       = 1000000000.0;   % = learning rate has 'blown up'
DEFAULT_BLOWUP_FAC   = 0.8;       % when lrate 'blows up,' anneal by this fac
block      = min(floor(5*log(frames)),floor(0.3*frames)); % heuristic 
lrate      = 0.00065/log(chans); 
annealdeg  = 60;
annealstep = .9;                      % defaults declared below
nochange   = 0.0000001;
maxsteps   = 512;

if ~exist('sphering','var')
    sphering = 'on';
end

if strcmp(sphering,'on'), %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
      weights = eye(chans); % begin with the identity matrix
  data = (sphere*data);      % actually decorrelate the electrode signals
elseif strcmp(sphering,'off') %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      sphere = 2.0*inv(sqrtm(cov(data'))); % find the "sphering" matrix = spher()
      weights = eye(chans)*sphere; % begin with the identity matrix
      sphere = eye(chans);                 % return the identity matrix
elseif strcmp(sphering,'none')
      weights = eye(chans); % begin with the identity matrix
  sphere = eye(chans,chans);
end
%
%%%%%%%%%%%%%%%%%%%%%%%% Initialize ICA training %%%%%%%%%%%%%%%%%%%%%%%%%
%
  lastt=fix((frames/block-1)*block+1);
  BI=block*eye(chans);
  degconst = 180./pi;
  startweights = weights;
  oldweights = startweights;
  lrates = zeros(1,maxsteps);
%
%%%%%%%% ICA training loop using the logistic sigmoid %%%%%%%%%%%%%%%%%%%
%
  step=0;
  blockno = 1;  % running block counter for kurtosis interrupts

  while step < maxsteps, %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      permute=randperm(frames); % shuffle data order at each step

    for t=1:block:lastt, %%%%%%%%% ICA Training Block %%%%%%%%%%%%%%%%%%%                                                      
       u=weights*data(:,permute(t:t+block-1)); 
       %%%%%%%%%%%%%%%%%%% Logistic ICA weight update %%%%%%%%%%%%%%%%%%%
       %y=1./(1+exp(-u));                                                %
       %weights=weights+lrate*(BI+(1-2*y)*u')*weights;                   %
       y = sign(u).*(1-exp(-abs(u)))./(1+exp(-abs(u)));
       dweight = lrate*(BI-(y*u'))*weights;
       %dweight = bsxfun(@times,dweight,exp(-1i*(diag(dweight))));
       weights = weights + dweight;
      blockno = blockno + 1;
    end % training block %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~mod(step,10)
        imagesc(complexIm(weights));drawnow;
    end
      oldwtchange = weights-oldweights;
      step=step+1; 
      lrates(1,step) = lrate;
      angledelta=0.;
      delta=reshape(oldwtchange,1,chans*chans);
      change=delta*delta'; 
      %%%%%%%%%%%%%%%%%%%% Anneal learning rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     if step> 2 
         angledelta=acos((delta*olddelta')/sqrt(change*oldchange));
              fprintf(...
         'step %d - lrate %5f, wchange %7.6f, angledelta %4.1f deg\n', ...
          step,lrate,change,degconst*angledelta);

     end
      oldweights = weights;
      if degconst*angledelta > annealdeg,  
        lrate = lrate*annealstep;          % anneal learning rate
        olddelta   = delta;                % accumulate angledelta until
        oldchange  = change;               %  annealdeg is reached
      elseif step == 1                     % on first step only
        olddelta   = delta;                % initialize 
        oldchange  = change;               
      end
%
%%%%%%%%%%%%%%%%%%%% Apply stopping rule %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
      if step >2 && change < nochange,      % apply stopping rule
           step=maxsteps;                  % stop when weights stabilize
      elseif change > DEFAULT_BLOWUP,      % if weights blow up,
        lrate=lrate*DEFAULT_BLOWUP_FAC;    % keep trying 
      end;                                 % with a smaller learning rate
      3;
end % end training %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  activations = weights*data;