function x  =getGradShift(data)

[u,s,v] = svds(data,1);
recon = u*s*v';
res_angles = circ_dist(angle(data),angle(recon));
xs = circ_std(res_angles,abs(data),[],2);
xm = circ_mean(res_angles,abs(data),2);
xm1 = angle(u);xm1 = circ_dist(xm1,circ_mean(xm1));
%x0 = mean(abs(res_angles));

%x0 = [circ_mean(angle(recon),abs(recon)); mean(abs(res_angles))/mean(abs(xm))];
%x = zeros(100,1);
%options = optimset('display','off');
x = res_angles'/xm1';
recon1 = recon.*exp(1i*xm1*x');
corrcoef(recon(9,:),data(9,:))
corrcoef(recon1(9,:),data(9,:))
figure;plot(abs(x));hold all;plot(circ_std(res_angles));
x  = filtfilt(gausswin(10),sum(gausswin(10)),x);
[u,s,v] = svds(data,2);
recon = u*s*v';
corrcoef(data(9,:),recon(9,:))
plot(abs(v(:,2))*200)
%res_angles = ci)rc_dist(angle(data),angle(recon1));
%x = res_angles'/xm1';
%plot(x);
%for i = 1:1000
%%    x(i,:) = lsqnonlin(@getCircDist,x0(:,i),[],[],[],data(:,i),xm1);
%x0 = res_angles(:,i)'/xm1';
%x(i) = lsqnonlin(@getRelDist,x0,[],[],options,res_angles(:,i),xm1);%,
%%scatter(xm1,res_angles(:,i),'filled');hold all;scatter(xm1,x(i)*xm1);hold off;drawnow;
%end
%scatter(x0(1:100),x)
%[u1 s1 v1] = svds(angle(data./recon),1);%-1%circ_dist(angle(data),angle(recon))
%xm2 = u1;

function err = getRelDist(x,data,bias)
err = circ_dist(data,bias*x);
%scatter(bias,data,'filled');hold all;scatter(bias,x*bias);hold off;drawnow;

function err = getCircDist(x,data,bias)
err = circ_dist(angle(data),x(2)*(bias-x(1)));
scatter(bias,angle(data));hold all;scatter(bias,angle(x(2)*bias-x(1)));hold off;drawnow;