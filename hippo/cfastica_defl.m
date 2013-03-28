function W = cfastica_defl(x)
eps = 0.1; % epsilon in G
[x,whiteningMatrix,dewhiteningMatrix,zerophaseMatrix] = whiten(x,.1);
n = size(x,1);

W = zeros(n);
maxcounter = 40;
for k = 1:n
    w = rand(n,1) + 1i*rand(n,1);
    counter = 0;
    wold = zeros(n,1);
    while min(sum(abs(abs(wold) - abs(w))), maxcounter - counter) > 0.001;
        wold = w;
        g = 1./(eps + abs(w'*x).^2);
        dg = -1./(eps + abs(w'*x).^2).^2;
        w = mean(x .* (ones(n,1)*conj(w'*x)) .* (ones(n,1)*g), 2) - ...
            mean(g + abs(w'*x).^2 .* dg) * w;
        w = w / norm(w);
        % Decorrelation:
        w = w - W*W'*w;
        w = w / norm(w);
        counter = counter + 1;
    end
    W(:,k) = w;
end
W = W';
x = W * x;
figure;imagesc(abs(x));