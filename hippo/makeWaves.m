function m = makeWaves(dims)
sigma = .3;
factor = 10;
parallel = 0;
linear = 0;
res = 10;
seed = 2*pi*rand(2,1);
envelope = ones(dims(3),2);
envelope = rand(dims(3),2);
envelope = filtfilt(gausswin(30),sum(gausswin(30)),envelope);
envelope = bsxfun(@minus,envelope,min(envelope));
envelope = bsxfun(@rdivide,envelope,max(envelope));
grid = repmat(linspace(0,2*pi,dims(1)),[dims(2) 1]);

if linear
    sFactor = factor;
else
    sFactor = 1;
end

if parallel
    grid1 = repmat(linspace(0,2*pi,dims(1))*sFactor,[dims(2) 1]);
else
    grid1 = repmat(linspace(0,2*pi,dims(1))'*sFactor,[1 dims(2)]);
end

m = sigma*randn([size(grid) dims(3)]);
for i = 1:dims(3)
    m(:,:,i) = m(:,:,i) + envelope(i,1)*exp(1i*(grid + i/res + seed(1)))+envelope(i,2)*exp(1i*(grid1+factor*i/res+seed(2)));
end