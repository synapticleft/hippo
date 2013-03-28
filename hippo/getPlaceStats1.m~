function getPlaceStats1(fields,thresh)
warning off;
fields(:) = zscore(fields(:));
fields = permute(fields,[1 3 2]);
buffer = 40;
%fields = cat(3,zeros(size(fields,1),size(fields,2),buffer),fields,zeros(size(fields,1),size(fields,2),buffer));
%allMax = [];
%counter = 0;
temp = squeeze(abs(mean(fields,3)));
threshInds = find(max(temp,[],2) > thresh);
fields = fields(threshInds,:,:);
lambda = 0;
betas = 0;
allLocs = [];
for i = 1:size(fields,1)
        filts{i} = [];
        tempFull = squeeze(fields(i,:,:));
        %tempFull = [zeros(size(tempFull,1),buffer) tempFull zeros(size(tempFull,1),buffer)];
        tempFullF = filtfilt(gausswin(12),sum(gausswin(12)),tempFull);
        tempF = mean(abs(tempFullF.'));
        [~,locs] = findpeaks(tempF,'minpeakheight',thresh,'minpeakdistance',buffer/2);
        allLocs =[allLocs locs];
        for j = 1:numel(locs)
            inds = max(1,locs(j)-buffer/2):min(size(fields,2),locs(j)+buffer/2);
            X = tempFull;
            X(~ismember(1:size(tempFull,1),inds),:) = 0;
            X = toeplitz(X(:),zeros(buffer/2,1));
            Y = fields;
            Y(i,inds,:) =0;
            %imagesc(abs(squeeze(Y(i,:,:))));return
            if locs(j) < size(fields,2)/2
                Y(:,end/2+1:end,:) = 0;
            else
                Y(:,1:end/2,:) = 0;
            end
            beta = (Y(:,:)*conj(X))/(X'*X + lambda*eye(size(X,2)));
            beta(beta < .1) =0;
            betas = betas+abs(beta);
            imagesc(abs(beta));title(i);pause(.1);
%             %temp = circshift(tempFull(:,locs(j)+(-buffer/2:buffer/2)),[0 0]);
%             winCorrs = zeros(buffer+1,size(fields,1));
%             for k = -buffer/2:buffer/2
%                 temp1 = fields(:,:,locs(j)+k+(-buffer/2:buffer/2));%tempFull(:,locs(j)+k+(-buffer/2:buffer/2));
%                 if locs(j) < size(fields,3)/2
%                     temp1(:,:,end/2+1:end) = 0;
%                 else
%                     temp1(:,:,1:end/2) = 0;
%                 end
%                 winCorrs(k+buffer/2+1,:) = (mean(bsxfun(@times,(temp1(:,:)),(temp(:)')),2));
%                 %filts{i}(k+buffer/2+1,j) = sum(temp(:).*conj(temp1(:)));
%             end
%             winCorrs(:,i) = 0;
%             [m1,m2] = max(abs(winCorrs));
%             [~,m3] = max(m1.*(m2 > buffer/2+1));
%             [~,m4] = max(m1.*(m2 < buffer/2+1));
%             subplot(121);hold off;imagesc(complexIm(winCorrs'));
%             if m1(m3) < thresh/2
%                 m3 = i;
%             end
%             if m1(m4) < thresh/2
%                 m4 = i;
%             end
%             hold all;scatter(m2(m4),m4,'filled');scatter(m2(m3),m3,'filled');%pause(1);
%             subplot(122);imagesc(imfilter(abs([squeeze(fields(i,:,:))' squeeze(fields(m3,:,:))' squeeze(fields(m4,:,:))']'),...
%                 fspecial('gaussian',5,1)));input('');%pause(1);
%             allMax{i}(:,j) =  [m2(m3) m2(m4) m1(m3) m1(m4)];
%             counter = counter + 1;
        end
end
%allMax = [allMax{:}];
%figure;scatter(allMax(1,:),allMax(3,:));hold all;scatter(allMax(2,:),allMax(4,:))
%f = [filts{:}];
%f = bsxfun(@rdivide,f,sqrt(sum(abs(f).^2)));