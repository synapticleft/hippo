function [y,yHat,yHats,acts,actsHat,yHats1,sts ] = regressSparse(stim,act,trials)
% Fit sparse activations using stimulus history

endAlign =0;
ord = 15;

X = [];y = [];Y = [];
for i = min(trials):max(trials)
    X = [X toeplitz(zeros(1,ord),(stim(trials == i)))];%cumsum
    Y = [Y toeplitz(zeros(1,ord),act(trials == i))];
    y = [y act(trials == i)];
end 
%X = zscore(X,0,2);
%y = [max(0,y);min(0,y)];
%y = zscore(y);
w = y/X;
yHat = w*X;

h = hist(trials,1:max(trials));
f = find(h);
sts = zeros(numel(f),max(h));acts = zeros(numel(f),max(h));actsHat = acts;elapse = [];


for i= 1:numel(f)
    if ~endAlign
        sts(i,1:h(f(i))) = stim(trials == f(i));
        acts(i,1:h(f(i))) = act(1,trials == f(i));
        actsHat(i,1:h(f(i))) = yHat(trials == f(i));
        elapse = [elapse 1:sum(trials == f(i))];
    else
        sts(i,end-h(f(i))+1:end) = stim(trials == f(i));
        acts(i,end-h(f(i))+1:end) = act(1,trials == f(i));
        actsHat(i,end-h(f(i))+1:end) = yHat(trials == f(i));
        elapse = [elapse sum(trials == f(i)):-1:1];
    end
    
end

% for i = 1:max(elapse)
%     XCov = X(:,elapse == i)*X(:,elapse == i)';
%     ws(i,:) = y(elapse == i)*X(:,elapse == i)'/(XCov + 1*diag(diag(XCov)));%y(elapse == i)/X(:,elapse == i);
%     yHats(elapse == i) = ws(i,:)*X(:,elapse == i);
%     %fs = find(elapse == i);
%     %yHats1(trials(fs),i) = yHats(elapse == i);
% end
% yHats1 = zeros(size(actsHat));
% for i= 1:numel(f)
%     if ~endAlign
%     yHats1(i,1:h(f(i))) = yHats(trials == f(i));
%     else
%         yHats1(i,end-h(f(i))+1:end) = yHats(trials == f(i));
%     end
% end

figure;plot(w')
figure;plot(y');hold all;plot(yHat');
corrcoef(y',yHat')