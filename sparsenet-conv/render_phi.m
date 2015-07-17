function array = render_phi(phi, m, dewhiteningMatrix)

[N J R] = size(phi);
% N1 = size(dewhiteningMatrix,1);
% if N1 > N
%     ratio = N1/N;
%     dewhiteningMatrix = dewhiteningMatrix((ratio-1)*N+(1:N),:);
% end
buf=1;

n = ceil(J/m);

array = -ones(buf+m*(N+buf),buf+n*(R+buf));

k = 1;

for i = 1:m
    for j = 1:n
        %bf = flipud(squeeze(phi(:,k,:)));
        bf = squeeze(phi(:,k,:));%dewhiteningMatrix*
        clim = max(abs(bf(:)));

        array(buf+(i-1)*(N+buf)+[1:N],buf+(j-1)*(R+buf)+[1:R]) = bf/clim;

        k = k+1;
        if k > size(phi,2)
            return
        end
    end
end

