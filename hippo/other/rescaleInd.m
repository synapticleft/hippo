function vec = rescaleInd(vec)

vec = zscore(vec);
margin = .05;
[vecS,inds] = sort(vec);

figure;plot(vecS);
i = input('lower bound (0 for none): ');
if i
    vec(inds(1:i)) = vecS(i+1)*(1+margin);%vec(inds(i+1))
end
i = input('upper bound (0 for none): ');

if i
    vec(inds(i:end)) = vecS(i-1)*(1+margin);%vec(inds(i-1))
end