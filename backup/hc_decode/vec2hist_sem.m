function [edges,ym,yl,yu] = vec2hist_sem(x,y,N,edges)

if nargin<4,
    edges = linspace(min(x)-10e-3,max(x)+10e-3,N);
end

[n,bin] = histc(x,edges);
for i=1:length(n)
    if sum(bin==i)>1
        ym(i) = mean(y(bin==i));
        yl(i) = std(y(bin==i))/sqrt(sum(bin==i));
        yu(i) = std(y(bin==i))/sqrt(sum(bin==i));
%         yl(i) = std(y(bin==i))/sqrt(sum(bin==i));
%         yu(i) = std(y(bin==i))/sqrt(sum(bin==i));
    else
        ym(i) = NaN; yl(i) = NaN; yu(i) = NaN;
    end
end