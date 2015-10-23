function dataOut = getKeyStrokes(file1,file2)

data = loadjson(file1);
data = data.data;

for i = 3:length(data)-1
    if ~iscell(data{i})
        dataOut(i-2,:) = data{i};
    else
        for j = 1:5
            if iscell(data{i}{j})
                if data{i}{j}{1}(1) == 'R'
                    dataOut(i-2,j) = -1*str2num(data{i}{j}{1}(2:end));
                elseif data{i}{j}{1} == 'L'
                    dataOut(i-2,j) = 0;
                else
                    dataOut(i-2,j) = str2num(data{i}{j}{1});
                end
            else
                dataOut(i-2,j) = data{i}{j};
            end
        end
    end
end
dataOut(dataOut(:,3) == 0,3) = 8;

% fid = fopen(['up_time' file2]);
% temp = textscan(fid,'%s %f-%f-%fT%f:%f:%f+01:00');
% bonsUp = getTime(temp);
% upKeys = temp{1};
% fclose(fid);

fid = fopen(['down_time' file2]);
temp = textscan(fid,'%s %f-%f-%fT%f:%f:%f+01:00');
bonsDown = getTime(temp);
downKeys = temp{1};
fclose(fid);
bonsDown = getKey(downKeys,bonsDown,dataOut(:,2));
figure;plot(bonsDown(2:end)-dataOut(dataOut(:,2) < 0,1))
dataOut(:,1) = dataOut(:,1) -bonsDown(1);

function out = getTime(struc)
out = struc{5}*3600+struc{6}*60+struc{7};

function bt = getKey(bk,bt,gk)
gk(gk >= 0) = [];
gk = gk*-1;
for i = 1:length(bk)
    if streq('D1',bk(i))
        break
    end
end
for j = 1:length(gk)
    if ~((gk(j) == 37 && streq(bk(i+j),'Left')) || (gk(j) == 39 && streq(bk(i+j),'Right')) ...
            || (gk(j) == 38 && streq(bk(i+j),'Up')) || (gk(j) == 40 && streq(bk(i+j),'Down')))
        'Error!!! key identities from bonsai and game dont match!!!'
        return
    end
end
bt = bt - bt(i);
bt(2:i) = [];
bt(j+2:end) = [];