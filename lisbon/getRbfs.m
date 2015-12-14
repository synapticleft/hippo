function [b,Xrbf ] = getRbfs(btype,n_basis,n_grid,X,s)
%GETRBFS Summary of this function goes here
%   Detailed explanation goes here

b.n=n;
b.type=btype;
b.mxv = linspace(0,1,n_basis(1));
b.myv = linspace(0,1,n_basis(2));
[xm,ym] = meshgrid(b.mxv,b.myv);
b.mxy = [xm(:) ym(:)];
if nargin>3
    b.s=s;
else
    b.s = diag([mean(diff(b.mxv)) mean(diff(b.myv))]).^2;
end

[xm,ym]=meshgrid(linspace(0,1,n_grid(1)),linspace(0,1,n_grid(2)));
b.x0=[xm(:) ym(:)];
b.rbf_basis = get2Dbasis_gaussian(b.x0,b.mxy,b.s);
if nargin>2
    Xrbf = get2Dbasis_gaussian(X,b.mxy,b.s);
end


function Xrbf = get2Dbasis_gaussian(X,rbfm,rbfs)

Xrbf=[];
for i=1:size(rbfm,1)
    Xrbf(:,i) = exp(-sum(bsxfun(@times,bsxfun(@minus,X,rbfm(i,:))*inv(rbfs),bsxfun(@minus,X,rbfm(i,:))),2)/2)/sqrt(det(rbfs))/(2*pi);
end

