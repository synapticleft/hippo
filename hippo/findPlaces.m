function [fits cis acts] = findPlaces(fields,thresh,sign)
% fit wavelets to place fields of complex-valued features

buffer = floor(size(fields,3)/2/10);
if ~exist('sign','var') || sign > 0
    sign = 1;
end
%[xs ys] = meshgrid(1:size(fields,2),1:size(fields,3));
for i = 1:size(fields,1)
    fits{i} = [];cis{i} = [];
    for j = 1:2
        tempFull = sign*squeeze(fields(i,:,(j-1)*size(fields,3)/2+(1:size(fields,3)/2)));
        %tempFull = filtfilt(gausswin(buffer/2),sum(gausswin(buffer/2)),tempFull')';
        temp = mean(tempFull);
        tempF = mean(abs(tempFull));
        [~,locs] = findpeaks(mean(tempF,1),'minpeakheight',thresh,'minpeakdistance',buffer);
        x0 = [];lb = [];ub = [];
        for k = 1:numel(locs)
            x0 = [x0; 1/20 -.1 locs(k) real(temp(locs(k)))' imag(temp(locs(k)))'];
            lb = [lb; 0 -inf locs(k)-buffer/2 -inf -inf];
            ub = [ub; inf 0 locs(k)+buffer/2 inf inf];
            %allMags(k,:) = tempF(:,locs(k))'/norm;
        end
        subplot(311);
        imagesc(complexIm(tempFull,[],[],[],[],thresh*2));title((2*(i-1)+j));
        hold all;scatter(locs,ones(numel(locs),1)*size(tempFull,2)/2,'w','filled');hold off;
        if numel(locs)
            %f = @(x,xdata) myfun(x,xdata,allMags);
            [beta,~,resid,~,~,~,J] = lsqcurvefit(@myfun,x0,size(temp,2),[real(temp);imag(temp)],lb,ub);
            ci = nlparci(beta,resid,'jacobian',J);
            ci = reshape(ci(:,2)-ci(:,1),size(beta));
            for k = 1:numel(locs)
                mor = makeCmor(size(temp,2),beta(k,1),beta(k,2),beta(k,3));
                acts{i}(k,:) = (mor*tempFull')/(mor*mor');
            end
            c = myfun(beta,size(temp,2),acts{i});%,allMags
            c = complex(c(1:size(c,1)/2,:),c(size(c,1)/2+1:end,:));
            subplot(312);imagesc(complexIm(c,[],[],[],[],thresh*2));hold all;
            scatter(beta(:,3),ones(size(beta,1),1)*size(tempFull,2)/2,'w','filled');hold off;title(numel(locs));drawnow;
            beta(:,3) = beta(:,3) + (j-1)*size(fields,3)/2;
            fits{i} = [fits{i} beta'];
            cis{i} = [cis{i} ci'];
            %fits(2*(i-1)+j).x
            tc = tempFull-c;
            tt = [mean(tempFull);mean(c)];
            subplot(313);sPlot(complex(abs(tt)/max(abs(tt(:)))*6,angle(tt)),[],0);%imagesc(complexIm(tc,[],[],[],[],thresh*2));
            drawnow;input('');%pause(1);
            %[sum(mean(tempFull.*conj(tempFull))) sum(mean(c.*conj(c))) sum(mean(tc.*conj(tc)))]
        end
        %ac = 0;
        %for k = 1:size(temp,1)
        %    ac = ac + xcov(temp(k,:),buffer);
        %end
        %
    end
end

function F = myfun(x,len,mags)%
F = 0;
for i = 1:size(x,1);
    thisMor = makeCmor(len,x(i,1),x(i,2),x(i,3));
    if exist('mags','var')
        F = F + mags(i,:)'*thisMor;%
    else
        F = F + complex(x(i,4),x(i,5))*thisMor;%
    end
end
F = [real(F);imag(F)];

function [psi,X] = makeCmor(len,Fb,Fc,off)
X = (1:len)-off;  % wavelet support.
psi = exp(2*1i*pi*Fc*X).*exp(-(X.*X)*Fb);%((pi*Fb)^(-0.5))*