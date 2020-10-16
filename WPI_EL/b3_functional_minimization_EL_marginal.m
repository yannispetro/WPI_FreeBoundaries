function [] = b3_functional_minimization_EL_marginal(... 
    MC_file, points, ord, targetPDFs, par, ti, tfs, Fss)
    
    ndof = par.ndof;
    options = bvpset('Vectorized','off','NMax',500);
    for jj = 1:length(tfs)
        tf = tfs(jj);
        Fs = Fss(jj);
        
        domain = c1_Get_bounds(MC_file,targetPDFs,tf,1e-4);
        [ODEs,Lagrangian,BC_free_sym] = EL_eqs(par,tf);
        [ tt,dt,N ] = Dim_time( Fs,tf );
        
        fprintf('%d | tf = %.5f | domain = ( %.5f,%.5f) \n', jj,tf,domain(1),domain(2))
        for kk = 1:length(targetPDFs)
            tic
            
            j = targetPDFs(kk);

            % Create the original grid
            X = linspace(domain(2*kk-1), domain(2*kk), points);
            BCs = Create_BCs(ndof,j,BC_free_sym);

            action = zeros([1, points]);
        %     parfor_progress(points);
            solinit = bvpinit(tt, zeros(2*ord*ndof, 1));
%             for r = 1:points
            for r = points:-1:1
                BVXj = X(r);

%                 solinit = bvpinit(tt, zeros(2*ord*ndof, 1));
                % Solve nonlinear BVP
%                 if jj == 1
%                     solinit = bvpinit(tt, zeros(2*ord*ndof, 1));
%                 else
%                     solinit = bvpinit(SOLS{kk,r}, [tt(1), tt(end)]);
%                 end
                sol = bvp4c(@M_ode, @M_bc, solinit, options, ODEs, BCs, BVXj);
%                 SOLS{kk,r} = sol;
                solinit = bvpinit(sol, [tt(1), tt(end)]);

                XX = deval(sol, tt);

                % Calculate action
                Integ = zeros(1,N);
                for it = 1:N
                    Integ(it) = Lagrangian([XX(:,it);tt(it)]);
                end

                action(r) = trapz(tt, Integ);
        %         parfor_progress;
            end
        %     parfor_progress(0);
            elapsed_time = toc;
            PDF = exp(-action);
            save(['files/WPI_EL_PDF' num2str(j) '_tf' num2str(tf) '_fs' num2str(Fs) '.mat'],...
            'points','domain','PDF','j','tf','par','Fs','elapsed_time');
        end

    end
end