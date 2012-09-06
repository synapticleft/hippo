function Phi = sparsenet(Xf,dw,zp)
% sparsenet.m - Olshausen & Field sparse coding algorithm
% 
% Before running you must first load the training data array IMAGES

num_trials=10000;
batch_size=1000;

% number of outputs
M= 30;

% number of inputs
N= size(Xf,1);

% initialize basis functions (comment out these lines if you wish to 
% pick up where you left off)
Phi=randn(N,M);
Phi=Phi*diag(1./sqrt(sum(Phi.*Phi)));

% learning rate (start out large, then lower as solution converges)
eta = .0005;

% lambda
lambda = .1;

a_var=ones(M,1);
var_eta=.1;

display_every=1;
Phi1 = zp*dw*Phi;
display_Phi(complex(Phi1(1:32,:),Phi1(34:end-1,:)),a_var,1);

for t=1:num_trials
    
    % choose an image for this batch

    inds=randperm(size(Xf,2));
    inds = inds(1:batch_size);
    I = Xf;%(:,inds);
    % calculate coefficients for these data via LCA

    ahat = sparsify(I,Phi,lambda,'hard');
    % calculate residual error

    R=I-Phi*ahat;

    % update bases

    dPhi = eta * R * ahat'; %  learning rule here
    Phi = Phi + dPhi;
    Phi=Phi*diag(1./sqrt(sum(Phi.*Phi))); % normalize bases

    % accumulate activity statistics
    a_var=(1-var_eta)*a_var + var_eta*mean(ahat.^2,2);
    % display
    if (mod(t,display_every)==0)
        Phi1 = zp*dw*Phi;
        display_Phi(complex(Phi1(1:32,:),Phi1(34:end-1,:)),a_var,1);
        sfigure(1);subplot(2,3,3);sPlot(ahat,[],0);
        subplot(2,3,6);bar(a_var);drawnow;
        %sfigure(1);subplot(2,3,3);imagesc(I);subplot(2,3,6);imagesc(I-R);drawnow;
    end
    
end

function array = display_Phi(A,Z_var,fig_num)

[L M]=size(A);
sz = [8 4];
if L == 96
    sz = [16 6];
end
%sz=sqrt(L);

buf=1;
m = ceil(sqrt(M));
n = ceil(M/m);

array=-ones(buf+n*(sz(1)+buf),buf+m*(sz(2)+buf))*(1+1j);
k=1;
for c=1:m
  for r=1:n
      if k > M
          break;
      end
    clim=max(abs(A(:,k)));
    array(buf+(r-1)*(sz(1)+buf)+[1:sz(1)],buf+(c-1)*(sz(2)+buf)+[1:sz(2)])=...
	reshape(A(:,k),sz(1),sz(2))/clim;
    k=k+1;
  end
end

subp_space = 0.03;
sfigure(fig_num);
clf;
colormap gray
subp(2,3,1,subp_space);
h(1)=imagesc(real(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('real')
subp(2,3,2,subp_space);
h(2)=imagesc(imag(array),[-1 1]);
axis image off
colormap gray; freezeColors
title('imag')
subp(2,3,4,subp_space);
h(3)=imagesc(abs(array),[0 max(abs(array(:)))]);
axis image off
title('abs')
colormap gray; freezeColors
subp(2,3,5,subp_space);
h(4)=imagesc(angle(array),[-pi pi]);
alpha(abs(array)/max(abs(array(:))));
axis image off
colormap hsv; freezeColors
title('angle')
drawnow;
% if exist('Z_var','var')  && ~isempty(Z_var) && (max(Z_var(:)) > 0)
%   sfigure(fig_num+2);
%   subplot(211)
%   bar(double(Z_var)), axis([0 N+1 0 double(max(Z_var))])
%   title('Z variance')
%   subplot(212)
%   normA=sqrt(sum(abs(A).^2));
%   bar(double(normA)), axis([0 N+1 0 double(max(normA))])
%   title('basis norm (L2)')
% end
