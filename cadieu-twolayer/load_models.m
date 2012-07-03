function load_models(sub)

curDir = ['/media/work/hippocampus/state/' sub];
f = dir([curDir '2*.*']);
Aold = zeros(96,20);
% for i = 1:numel(f)
%     if i > 1
%         Aold = m.A;
%     end
%     load([curDir f(i).name],'m');
%     m.A = bsxfun(@times,m.A,exp(1i*-angle(mean(m.A))));
%     display_A(m);pause(.01);
%     if i > 1
%         figure(99);scatter(i,log10(1-mean(abs(sum(m.A.*conj(Aold))))),'filled','b');hold on;
%         scatter(i,log10(1-min(abs(sum(m.A.*conj(Aold))))),'filled','r');hold on;
%     end
% end
counter = 0;
AAll = zeros(96,64*20);
for i = 1:numel(f)
    load([curDir f(i).name],'m');
    if i == numel(f) || ~strcmp(f(i).name(1:20),f(i+1).name(1:20))
        AAll(:,counter*20+(1:20)) = bsxfun(@times,m.A,exp(1i*-angle(mean(m.A))));
        counter = counter+1;
        %display_Ahelper(bsxfun(@times,m.A,exp(1i*-angle(mean(m.A)))),counter);counter = counter + 1;
    end
end
[a b c] = ica_lfp_complex(AAll,20);
complexMovies(a,[],[16 6]);
