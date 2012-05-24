function [x, V, A,R] = complexKalman(spikes,zSpikeFit,coeff,noise,trial,ref)

spikes = spikes(trial,:);
ref = ref(trial,:);
p = numel(coeff)-1;
ref(2:p,:) = 0;
for i = 2:p
    ref(i,i:end) = ref(1,1:(end-i+1));
end
A = [-coeff(2:end);eye(p)];A(end,:) = [];
% [x,y] = eig(A);
% y = mean(diag(y));
% for i = 1:(p-1)
%     zSpikeFit.mu(i+1) = zSpikeFit.mu(i)/y;
% end
zSpikeFit.mu = zSpikeFit.mu;
E = zeros(p); E(1,1) = noise;
% P = zeros(p);
% for i = 1:numIt
%     P = A*P*A' + E;
% end
init_x = zeros(p,1);
init_V = getCov(coeff,noise,p);
warning off all;
%vecP = (eye(p^2,p^2) - kron(conj(A),A))\E(:);
%P = reshape(vecP,p,p);
%R = P*zSpikeFit.Sigma(1,1)/P(1,1);
R = zSpikeFit.Sigma;
y = zSpikeFit.mu*ones(size(spikes));
y(:,~spikes) = 0;

%y = getState(ref,p);
warning on all;
tic;
[x,V] = myReset(y,A,E,R,init_x,init_V,spikes);%OneKalman
toc

function r = makeRand(s,R)
r = (complex(randn(s),randn(s))'*sqrtm(R)).'/sqrt(2);

function z1 = getState(z,p)
z1 = zeros(p,numel(z));
for i = 1:p
    z1(i,i:end) = z(1:(end-i+1));
end

function c = getCov(coeff,noise,p)
thisNoise = makeComplex(10^6,sqrt(noise),[1 0]);
z = arData(thisNoise,coeff);
for i = 1:p
    temp(i,:) = z(i:(end-p+i));
end
temp = flipud(temp);
c = cov(temp.');