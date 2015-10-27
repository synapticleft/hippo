function [a,d] = adStructToMatrix(signal)
%% convert the acquisition board signal struct into 2 matrices for analog and digital signals

%names = fieldnames(signal);

for i = 0:10
    try
        a(i+1,:) = eval(['signal.A' num2str(i)]);
    catch
        break
    end
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

a = a(:,1:(size(d,2)*10));
a = double(a) - 2^15;
d = double(d);