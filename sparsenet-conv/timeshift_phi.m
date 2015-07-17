function phi = timeshift_phi(phi,totResp)

[N J R] = size(phi);

%center = 'max';%'mass';

%mid = ceil(R/2);

%new_phi = zeros(size(phi));
frac = round(R*.5);
for j = 1:J
%     switch center
%         case 'mass'
%             den = sum(sum(abs(phi(:,j,:))));
%             ind = ceil(  sum([1:R]'.*squeeze(sum(abs(phi(:,j,:)), 1))) / den );
%             delta = mid - ind;
%     if delta <= 0
%         new_phi(:,j,1:end+delta) = phi(:,j,1-delta:end);
%     elseif delta > 0
%         new_phi(:,j,1+delta:end) = phi(:,j,1:end-delta);
%     end
%         case 'max'
            phi(:,j,:) = sign(totResp(j))*phi(:,j,:);
            [~,mx] = max(squeeze(mean(abs(phi(:,j,:)),1)));
            %temp = squeeze(phi(:,j,:));
            %[~,mx] = find(abs(temp) == max(abs(temp(:))));
            %[~,mx] = max(abs(phi(1,j,:)));
            phi(:,j,:) = circshift(phi(:,j,:),[0 0 frac-mx]);    
            %    end
end