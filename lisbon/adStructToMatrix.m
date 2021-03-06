function [a,d] = adStructToMatrix(signal,interGSR)
%% convert the acquisition board signal struct into 2 matrices for analog and digital signals

%names = fieldnames(signal);
a = [];d = [];
for i = 0:10
    try
        a(i+1,:) = eval(['signal.A' num2str(i)]);
    catch
        break
    end
end

if exist('interGSR','var')
    a = a(:,1:(floor(size(a,2)/interGSR)*interGSR));
    a = reshape(a,size(a,1),interGSR,[]);
end

for i = 0:10
    try
        d(3*(i)+1,:) = eval(['signal.D' num2str(i) 'X']);
        d(3*(i)+2,:) = eval(['signal.D' num2str(i) 'Y']);
        d(3*(i)+3,:) = eval(['signal.D' num2str(i) 'Z']);
    catch
        break
    end
end

if numel(d)
    a = a(:,1:(size(d,2)*10));
end
a = double(a) - 2^15;

d = double(d);