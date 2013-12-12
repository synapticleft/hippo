
function pall = decodeBayesian_gauss(y,lam,sigma,od)

if nargin<4, od=0; end

pall=zeros(size(y,1),size(lam,1));
for i=1:size(y,1)
    for p=1:size(lam,1)
        P(:,p) = normpdf(y(i,:)',lam(p,:)',sigma)+od;
%         keyboard
%         P(:,p) = 1/sqrt(2*pi)./sigma.*exp(-(y(i,:)'-lam(p,:)').^2/2./sigma.^2);
    end
    P = bsxfun(@rdivide,P,sum(P,2));
    pall(i,:) = sum(log(P));
%     keyboard
%     imagesc(P)
%     drawnow
%     pause(0.001)
%     keyboard
end

pall=(exp(bsxfun(@minus,pall,max(pall))));
pall=bsxfun(@rdivide,pall,sum(pall,2));

pall(~isfinite(sum(pall,2)),:)=1/size(pall,2);