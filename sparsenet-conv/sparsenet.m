
for t = 1:num_trials

    %% select data
    
    j = ceil(max(reg)*rand);
    Xsamp = X(:,reg == j);
    %% compute the map estimate
    tic
    S = size(Xsamp,2);
    P = S+R-1;	% number of selection locations for source generation
    a0 = zeros(J, P);
    %% no bounds
    lb  = zeros(1,J*P); % lower bound
    ub  = zeros(1,J*P); % upper bound
    nb  = ones(1,J*P); % bound type (none)
    
    [a1,fx,exitflag,userdata] = lbfgs(@objfun_a_conv, a0(:), lb, ub, nb, opts_lbfgs_a, Xsamp, phi, lambda);
    a1 = reshape(a1, J, P);
    time_inf = toc;

    %% reconstruct
    EI = zeros(N,S);
    for r = 1:R
        EIr = phi(:,:,r)*a1;
        srt = R+1-r;
        fin = R+S-r;
        EI = EI + EIr(:,srt:fin);
    end
    totResp = totResp + sum(a1,2);
    %% compute snr
    E = Xsamp - EI;
    snr = 10 * log10 ( sum(Xsamp(:).^2) / sum(E(:).^2) );

    switch mintype_lrn
        case 'minimize'
            [obj0,g] = objfun_phi(phi(:), Xsamp, a1, gamma, EI, dphi);

            phi1 = minimize(phi(:), 'objfun_phi', lrn_searches, Xsamp, a1, gamma, EI, dphi);
            phi1 = reshape(phi1,N,J,R);

            [obj1,g] = objfun_phi(phi1(:), Xsamp, a1, gamma, EI, dphi);

            phi = phi1;

        case 'gd'

            [obj0,g] = objfun_phi(phi(:), Xsamp, a1);
            dphi = reshape(g, N, J, R);

            phi1 = phi - eta * dphi;

            [obj1,g] = objfun_phi(phi1(:), Xsamp, a1);
            
            % pursue a constant change in angle
            angle_phi = acos(phi1(:)' * phi(:) / sqrt(sum(phi1(:).^2)) / sqrt(sum(phi(:).^2)));
            if angle_phi < target_angle
                eta = eta*eta_up;
            else
                eta = eta*eta_down;
            end

            if obj1 > obj0
                fprintf('objective function increased\n');
            else
                phi = phi1;
            end
    end


    %% truncate the log
    eta_log = eta_log(1:update-1);
    %% append
    eta_log = [ eta_log ; eta ];

    %% renormalize basis functions to have unit length
    for j = 1:J
        phi(:,j,:) = phi(:,j,:) / sqrt(sum(sum(phi(:,j,:).^2)));
    end

    %% compute the objective function after renormalization
    [obj2,g] = objfun_phi(phi(:), Xsamp, a1);
    
    fprintf('%s up %06d', paramstr, update);
    fprintf(' obj0 %.4f obj1 %.4f obj2 %.4f ang %.4f', ...
            obj0, obj1, obj2, angle_phi);
    fprintf(' snr %.4f eta %.8f inf %.4f\n', snr, eta, time_inf);

    if display_every == 1 || mod(t,display_every) == 1

        %% display image, reconstruction, error
        figure(10); clf;
        
        mn = min([Xsamp(:) ; EI(:) ; E(:)]);
        mx = max([Xsamp(:) ; EI(:) ; E(:)]);
        subplot(3,1,1); imagesc(Xsamp, [mn mx]); %axis image off; colorbar;
        subplot(3,1,2); imagesc(EI, [mn mx]); %axis image off; colorbar;
        subplot(3,1,3); imagesc(E, [mn mx]); %axis image off; colorbar;

        %% display coefficients
        figure(6); if min(size(a1)) == 1 plot(a1); else imagesc(a1); end

        %% display basis functions
        array = render_phi(phi, Jrows,dewhiteningMatrix); %render_phi_2d(phi,Jrows);

        figure(7); imagesc(array);% axis image off; colormap(gray);

        %% plot our dynamic eta
        figure(8);
        plot(eta_log)

        drawnow;
    end

    
    if (update == 1 || t == num_trials)%(save_every == 1 || mod(update,save_every) == 1)

        %% make the output directory
        [sucess,msg,msgid] = mkdir(sprintf('state/%s', paramstr));

        %% write the basis functions as a png
        %% max(abs(array)) <= 1 and there are 64 entries in a colormap, so..
        imwrite((array+1)*32, colormap, ...
            sprintf('state/%s/bf_up=%06d.png',paramstr,update), ...
            'png');
    
        %% save the basis functions as a matlab variable
        eval(sprintf('save state/%s/phi.mat phi', paramstr));

        %% save all other useful parameters to a separate matlab file
        saveparamscmd = sprintf('save state/%s/params.mat', paramstr);
        saveparamscmd = sprintf('%s lambda', saveparamscmd);
       % saveparamscmd = sprintf('%s gamma', saveparamscmd);
        saveparamscmd = sprintf('%s dewhiteningMatrix', saveparamscmd);
        saveparamscmd = sprintf('%s whiteningMatrix', saveparamscmd);
        saveparamscmd = sprintf('%s eta', saveparamscmd);
        saveparamscmd = sprintf('%s eta_up', saveparamscmd);
        saveparamscmd = sprintf('%s eta_down', saveparamscmd);
        saveparamscmd = sprintf('%s eta_log', saveparamscmd);
        saveparamscmd = sprintf('%s J', saveparamscmd);
        saveparamscmd = sprintf('%s R', saveparamscmd);
       % saveparamscmd = sprintf('%s datatype', saveparamscmd);
        saveparamscmd = sprintf('%s mintype_inf', saveparamscmd);
        saveparamscmd = sprintf('%s mintype_lrn', saveparamscmd);
        saveparamscmd = sprintf('%s update', saveparamscmd);
       % saveparamscmd = sprintf('%s reload_every', saveparamscmd);
        eval(saveparamscmd);

        %% save the coefficient and image/reconstruction/error figures
        saveas(6, sprintf('state/%s/activation.png', paramstr));
        saveas(10, sprintf('state/%s/reconstruction.png', paramstr));
    end


    update = update + 1;
end