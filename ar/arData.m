function y = arData(x,b)
y = x;
if numel(x) == max(size(x))
    for i = numel(b):numel(x)
        for j = max(1,i-numel(b)+1):(i-1)
            y(i) = y(i) - y(j)*b(i-j+1);
        end
    end
else
    for k = 1:size(x,1)
        for i = numel(b):size(x,2)
            for j = max(1,i-numel(b)+1):(i-1)
                y(k,i) = y(k,i) - y(k,j)*b(i-j+1);
            end
        end
    end
end