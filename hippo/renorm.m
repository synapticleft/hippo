function x = renorm(x)

for i = 1:size(x,1)
    x(i,:) = x(i,:) - min(x(i,:));
    x(i,:) = x(i,:)/max(x(i,:));
end