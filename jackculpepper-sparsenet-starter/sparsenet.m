for t = 1:num_trials


    switch datasource
        case 'images'
            %% choose an image for this batch
            i = ceil(K*rand);
            I = IMAGES(:,:,i);

            X = zeros(L,B);

            % extract subimages at random from this image to make data vector X
            for b = 1:B
                r = buff + ceil((Nsz-Lsz-2*buff)*rand);
                c = buff + ceil((Nsz-Lsz-2*buff)*rand);

                X(:,b) = reshape(I(r:r+Lsz-1,c:c+Lsz-1), L, 1);

                %% normalize each patch to mean zero, unit variance
                %% note: you may not want to do this! however, doing so makes
                %% it easy to get things working on a variety of data sets.

                X(:,b) = X(:,b) - mean(X(:,b));
                X(:,b) = X(:,b) / std(X(:,b));
            end
        case 'lfp'
            temp = randperm(size(Xf,2));
            X = Xf(:,temp(1:B));
    end
    
    tic
    switch mintype_inf
        case 'lbfgsb'

            %% stupidest initialization, but consistent
            a0 = zeros(M, 1);
            a = zeros(M, B); 

            %% run lbfgsb on each patch, separately.

            %% note this could be sped up by running on the entire batch.
            %% however, for large batches you may then need to increase the
            %% number of l-bfgs-b iterations to get a consistent solution.

            %% running separately is slower, but consistent.

            for b = 1:B
                [a1,fx,exitflag,userdata] = lbfgs(@objfun_a,a0(:),lb,ub,nb, ...
                    opts,phi,X(:,b),lambda);
                a(:,b) = a1;
            end
        case 'batch'
            
            a = zeros(M,B);%phi\X;
            [a] = lbfgs(@objfun_a,a(:),lb,ub,nb,opts,phi,X,lambda);
            a = reshape(a,M,B);
    end
    time_inf = toc;
    E = X - phi*a;
    snr = 10 * log10 ( sum(X(:).^2) / sum(E(:).^2) );

    %snr_log = snr_log(1:update-1);
    snr_log = [ snr_log ; snr ];


    tic
    % update bases
    switch mintype_lrn
        case 'gd'
            [obj0,g] = objfun_phi(phi(:),a,X,lambda,gamma);
            dphi = reshape(g,L,M);

            phi1 = phi - eta*dphi;

            [obj1,dphi] = objfun_phi(phi1(:),a,X,lambda,gamma);

            if obj1 > obj0
                fprintf('warning: objfun increased\n');
            end

            %% pursue a constant change in angle
            angle_phi = acos(phi1(:)' * phi(:) / sqrt(sum(phi1(:).^2)) / sqrt(sum(phi(:).^2)));
            if angle_phi < target_angle
                eta = eta*1.01;
            else
                eta = eta*0.99;
            end

            phi = phi1;

            eta_log = eta_log(1:update-1,:);
            eta_log = [ eta_log ; [eta angle_phi]];


    end
    time_lrn = toc;


    if test_every == 1 || mod(update,test_every) == 0
        %% do inference on the test set
        switch mintype_inf
            case 'lbfgsb'
                atest1 = zeros(M, Btest);
                for i = 1:Btest
                    atest0 = zeros(M,1);
                    [atest1(:,i),fx,exitflag,userdata] = lbfgs(@objfun_a, ...
                        atest0,lb,ub,nb,opts,phi,Xtest(:,i),lambda);
                    fprintf('\r%d / %d', i, Btest);
                end
                %atest0 = atest1;
            case 'batch'
                atest1 = zeros(M,Btest);%phi\Xtest;
                %lb  = zeros(1,numel(atest1)); % lower bound
                %ub  = zeros(1,numel(atest1)); % upper bound
                %nb  = ones(1,numel(atest1));  % bound type (lower only)
                %nb  = zeros(1,numel(atest1)); % bound type (none)
                [atest1] = lbfgs(@objfun_a,atest1(:),lb,ub,nb,opts,phi,Xtest,lambda);
                atest1 = reshape(atest1,M,numel(atest1)/M);
        end
        fprintf('\n');
        objtest = objfun_a(atest1(:),phi,Xtest,lambda);
        objtest_log = [ objtest_log objtest ];

        sfigure(7);
        plot(test_every*(1:length(objtest_log)), objtest_log, 'r-');
        title('Test Set Energy History');
        xlabel('Iteration');
        ylabel('E');
    end


    %% display

    if display_every == 1 || mod(update,display_every) == 0
        % Display the bfs
        %array = render_network(phi, Mrows);
%           sfigure(1); colormap(gray);
%        imagesc(array, [-1 1]);
%        axis image off;

%        phi1 = zerophaseMatrix*dewhitenMatrix*phi;
%        display_Phi(complex(phi1(1:size(phi1,1)/2-1,:),phi1(size(phi1,1)/2+1:end-1,:)),1);

 
        EI = phi*a;

        mx = max(abs([ EI(:) ; X(:) ]));
        mn = min(abs([ EI(:); X(:)]));
        sfigure(4);
        subplot(2,2,1),imagesc(EI, [mn mx]),title('EI');
            colormap(gray),axis image off;
        subplot(2,2,2),imagesc(X,[mn mx]),title('X');
            colormap(gray),axis image off;
        subplot(223);%imagesc(abs(a)),title('a');colormap(gray),axis image off;
        a = bsxfun(@rdivide,abs(a),max(abs(a)));
        [h,x] = hist(a(:),0:.1:1);
        plot(x,min(2,h/size(a,2)));
        xs = meshgrid(1:B,1:M);
        subplot(224);
        scatter(xs(:),abs(a(:)),'filled');hold all;
        scatter(xs(1,:),max(abs(a)),'filled','r');hold off;
        sfigure(5);
        bar(std(a,0,2));
        axis tight;

        sfigure(6);
        plot(1:update, eta_log, 'r-');
        title('\eta History');
        xlabel('Iteration');
        ylabel('\eta');

        sfigure(12);
        plot(1:update, snr_log, 'r-');
        title('Reconstruction SNR History');
        xlabel('Iteration');
        ylabel('SNR');

        drawnow;

    end
        if save_every == 1 || mod(update,save_every) == 0
            phis(floor(update/save_every),:,:) = phi;
%             array_frame = uint8(255*((array+1)/2)+1);
% 
%             [sucess,msg,msgid] = mkdir(sprintf('state/%s', paramstr));
%  
%             imwrite(array_frame, ...
%                 sprintf('state/%s/phi_up=%06d.png',paramstr,update), 'png');
%             eval(sprintf('save state/%s/phi.mat phi',paramstr));
% 
%             saveparamscmd = sprintf('save state/%s/params.mat', paramstr);
%             saveparamscmd = sprintf('%s lambda', saveparamscmd);
%             saveparamscmd = sprintf('%s gamma', saveparamscmd);
%             saveparamscmd = sprintf('%s eta', saveparamscmd);
%             saveparamscmd = sprintf('%s eta_log', saveparamscmd);
%             saveparamscmd = sprintf('%s objtest_log', saveparamscmd);
%             saveparamscmd = sprintf('%s L', saveparamscmd);
%             saveparamscmd = sprintf('%s M', saveparamscmd);
%             saveparamscmd = sprintf('%s mintype_inf', saveparamscmd);
%             saveparamscmd = sprintf('%s update', saveparamscmd);
% 
%             eval(saveparamscmd);
    

        end

    %% renormalize columns of phi to have unit length
    phi = phi*diag(1./sqrt(sum(phi.^2)));

    fprintf('%s', paramstr);
    fprintf(' update %d o0 %.8f o1 %.8f eta %.4f', update, obj0, obj1, eta);
    fprintf(' ang %.4f', angle_phi);
    fprintf(' snr %.4f inf %.4f lrn %.4f\n', snr, time_inf, time_lrn);

    update = update + 1;
end
