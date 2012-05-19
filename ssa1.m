function [outSer U V]= ssa1(x1,L,sv)%[y,r,vr]=U1 [U V u] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -----------------------------------------------------------------                           
%    Author: Francisco Javier Alonso Sanchez    e-mail:fjas@unex.es
%    Departament of Electronics and Electromecanical Engineering
%    Industrial Engineering School
%    University of Extremadura
%    Badajoz
%    Spain
% -----------------------------------------------------------------
%
% SSA generates a trayectory matrix X from the original series x1
% by sliding a window of length L. The trayectory matrix is aproximated 
% using Singular Value Decomposition. The last step reconstructs
% the series from the aproximated trayectory matrix. The SSA applications
% include smoothing, filtering, and trend extraction.
% The algorithm used is described in detail in: Golyandina, N., Nekrutkin, 
% V., Zhigljavsky, A., 2001. Analisys of Time Series Structure - SSA and 
% Related Techniques. Chapman & Hall/CR.

% x1 Original time series (column vector form)
% L  Window length
% y  Reconstructed time series
% r  Residual time series r=x1-y
% vr Relative value of the norm of the approximated trajectory matrix with respect
%	  to the original trajectory matrix

% The program output is the Singular Spectrum of x1 (must be a column vector),
% using a window length L. You must choose the components be used to reconstruct 
%the series in the form [i1,i2:ik,...,iL], based on the Singular Spectrum appearance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pca = 1;
[u s v] = svds(x1,sv);
%figure;plot(log(diag(s)));
v = s*v';
v = x1;
% Step1 : Build trayectory matrix
   N=size(v,2);%length(x1); 
   if L>N/2;L=N-L;end
	K=N-L+1;
   X=zeros(L*size(v,1),K);  
	for i=1:K
	  temp=v(:,i:L+i-1); 
      X(:,i) = temp(:);
    end
    %X1 = toeplitz(toeplitz([x1(1,1); zeros(size(x1,1)*L-1,1)],x1);
    %X1 = flipud(X1(:,L:end));
    
% Step 2: SVD
if pca
   S=X*X';
	[U,d]=eigs(S,min(size(S,1)-1,100));
    d=diag(d);
	%[d,i]=sort(-diag(autoval));  
   %d=-d;
   %U=U(:,i);
   sev=sum(d); 
	figure;plot(log10(d/sev)),hold on,plot(log10(d/sev),'rx');
	title('Singular Spectrum');xlabel('Eigenvalue Number');ylabel('Eigenvalue (% Norm of trajectory matrix retained)')
   V=(X')*U;
   %rc=U*V';

% Step 3: Grouping

   I=input('Choose the components to reconstruct the series: ');
   if numel(I) == 1 I = 1:I; end
   %Vt=V';
   rca=U(:,I)*V(:,I)';%Vt(I,:);
       params.Fs = 1250/32;
    [S,f] = mtspectrumc(V,params);
    sPlot(sqrt(S)',f);
else
    [a,b,c,~,e,~,g] = runica(X);%,'extended',-2
    e
    figure;plot(c);
    params.Fs = 1250/32;
    [S,f] = mtspectrumc(g',params);
    sPlot(sqrt(S)',f);
    I=input('Choose the components to reconstruct the series: ');
    if numel(I) == 1 I = 1:I; end
    rca = b(:,I)\(a(I,:)\g(I,:));
    size(rca)
end

   
%    U1 = reshape(U,[sv L size(U,2)]);
% figure;
% for i= 1:size(U1,1)
%     subplot(1,size(U1,1),i);sPlot(squeeze(U1(i,:,:))',1:size(U1,2),0);
%     %imagesc(squeeze(U1(i,:,:)));
% end

% Step 4: Reconstruction

   y=zeros(size(v));%(N,1);  
   Lp=min(L,K);
   Kp=max(L,K);
for i = 1:size(y,1)
    temp = rca(i:size(y,1):end,:);
   for k=0:Lp-2
     for m=1:k+1;
      y(i,k+1)=y(i,k+1)+(1/(k+1))*temp(m,k-m+2);
     end
   end

   for k=Lp-1:Kp-1
     for m=1:Lp;
      y(i,k+1)=y(i,k+1)+(1/(Lp))*temp(m,k-m+2);
     end
   end

   for k=Kp:N
      for m=k-Kp+2:N-Kp+1;
       y(i,k+1)=y(i,k+1)+(1/(N-k))*temp(m,k-m+2);
     end
   end
end
outSer = u*s*y;
% 
%    figure;subplot(2,1,1);hold on;xlabel('Data poit');ylabel('Original and reconstructed series')
%    plot(x1);grid on;plot(y,'r')
%    
%    r=x1-y;
%    subplot(2,1,2);plot(r,'g');xlabel('Data poit');ylabel('Residual series');grid on
%    vr=(sum(d(I))/sev)*100;
