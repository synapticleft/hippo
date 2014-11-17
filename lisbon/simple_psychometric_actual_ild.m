



bad_index=find(responses==3);
responses(bad_index)=[];
times_in_sound(bad_index)=[];
ILD_timecourses(bad_index,:)=[];


    ILD=[-3 -2 -1 1 2 3];
    ZeroPadding =.1;
    num_of_chunks_max=stim_duration_sec/(chunk_duration_msec/1000);
    
 actual_timecourses=nan(size(ILD_timecourses));   
for i =1:length(times_in_sound)
 actual_timecourses(i,(1:min(floor((times_in_sound(i)-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max))) = ILD_timecourses(i,(1:min(floor((times_in_sound(i)-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)));
actual_means=nanmean(actual_timecourses,2);
end
    
  mean(ILD_timecourses(1,1:min(floor((times_in_sound(1)-.1-length(ZeroPadding)/Fs-ramp_stim_dur_msec/1000)/(chunk_duration_msec/1000)),num_of_chunks_max)))
[a b] = hist(actual_means,6);

center_ILDs=b;
b=b+mean(diff(b))/2;

for i = 1:length(times_in_sound)
    
    if actual_means(i)<=b(1)
        conditions_order(i)=1;
    elseif actual_means(i)>b(1) && actual_means(i)<=b(2)
        conditions_order(i)=2;
    elseif actual_means(i)>b(2) && actual_means(i)<center_ILD
        conditions_order(i)=3;
    elseif actual_means(i)>center_ILD && actual_means(i)<=b(4)
        conditions_order(i)=4;
    elseif actual_means(i)>b(4) && actual_means(i)<b(5)
        conditions_order(i)=5;
    elseif actual_means(i)>=b(5) 
        conditions_order(i)=6;
    end
end
    
cond_1=find(conditions_order==1);
one_1=find(responses(cond_1)==1);
p_1=length(one_1)/length(cond_1);

cond_2=find(conditions_order==2);
one_2=find(responses(cond_2)==1);
p_2=length(one_2)/length(cond_2);

cond_3=find(conditions_order==3);
one_3=find(responses(cond_3)==1);
p_3=length(one_3)/length(cond_3);

cond_4=find(conditions_order==4);
one_4=find(responses(cond_4)==1);
p_4=length(one_4)/length(cond_4);

cond_5=find(conditions_order==5);
one_5=find(responses(cond_5)==1);
p_5=length(one_5)/length(cond_5);

cond_6=find(conditions_order==6);
one_6=find(responses(cond_6)==1);
p_6=length(one_6)/length(cond_6);


p_all=[p_1 p_2 p_3 1-p_4 1-p_5 1-p_6];

p_all = [p_all; p_all];

p_mean=mean(p_all);


%g_1=[1;1;1;1;1;2;2;2;2;2];

%p=anovan(p_combined, g_1)

    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fit SIGMOIDAL extracted from cftool%%%%%%
ILD=center_ILDs;
 
 ok_ = isfinite(ILD) & isfinite(p_mean);
if ~all( ok_ )
    warning( 'GenerateMFile:IgnoringNansAndInfs',...
        'Ignoring NaNs and Infs in data.' );
end
 
starting_point = [.5 .5];   %extract coefficients from cftool. 
sigmoidal_fit = fittype('1 / (1 + exp(a*(b + x)))',...
    'dependent',{'y'},'independent',{'x'},...
    'coefficients',{'a', 'b'});

% Fit this model using new data
sigmoid_curve = fit(ILD(ok_)',p_mean(ok_)',sigmoidal_fit,'Startpoint',starting_point);


 E = std(p_all).*ones(size(ILD)); 

plot(sigmoid_curve, 'black');



hold on


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

plot(ILD(1),p_mean(1),'o', 'MarkerSize', 12, 'LineWidth', 2,'color',[1 0 0]);
plot(ILD(2),p_mean(2),'o', 'MarkerSize', 12,'LineWidth', 2,'color',[1 0.5 0]);
plot(ILD(3),p_mean(3),'o', 'MarkerSize', 12, 'LineWidth', 2,'color',[1 .75 0]);
plot(ILD(4),p_mean(4),'o', 'MarkerSize', 12,'LineWidth', 2,'color',[0 .75 1]);
plot(ILD(5),p_mean(5),'o', 'MarkerSize', 12, 'LineWidth', 2, 'color',[0 0.5 1]);
plot(ILD(6),p_mean(6),'o', 'MarkerSize', 12,'LineWidth', 2, 'color',[0 0 1]);

errorbar(ILD,p_mean,E,'.k', 'MarkerSize', .1)
xlim([ILD(1)-.1 ILD(6)+.1])
ylim([-.05 1.05])
  
  
  title('Behavioral Performance','FontSize',15)
  ylabel('Proportion Left','FontSize',15)
  xlabel('Average ILD (dB)','FontSize',15')
  
 set(gca,'YTick', [0 .25 .5  .75 1])
 
 set(gca,'YTickLabel','0|0.25|0.5|0.75|1')
  
  curve_values=sigmoid_curve(-6:.0001:6);
  ILD_range=-6:.0001:6;
  
  point_1=find(curve_values>.89999 & curve_values<.9001,1,'last');
  sILD_1=ILD_range(point_1);
  
   point_2=find(curve_values>.59999 & curve_values<.6001,1,'last');
   sILD_2=ILD_range(point_2);
   
    point_3=find(curve_values>.49999 & curve_values<.5001,1,'last');
    sILD_3=ILD_range(point_3);
    
     point_4=find(curve_values>.39999 & curve_values<.4001,1,'last');
     sILD_4=ILD_range(point_4);
     
      point_5=find(curve_values>.09999 & curve_values<.1001,1,'last');
      sILD_5=ILD_range(point_5);

