function [] = b3_functional_minimization_EL_joint(... 
    MC_file, points, ord, targetPDF, par, ti, tfs, Fss)
    
    ndof = par.ndof;

    options = bvpset('Vectorized','off');
    
    j1 = targetPDF(1);
    j2 = targetPDF(2);
    
    for jj = 1:length(tfs)
        tic
        tf = tfs(jj);
        Fs = Fss(jj);
    
        domain = c1_Get_bounds(MC_file,targetPDF,tf,1e-4);
        [ODEs,Lagrangian,BC_free_sym] = EL_eqs(par,tf);

        % time
        [ tt,dt,N ] = Dim_time( Fs,tf );



        % Create the original grid 
        x1 = linspace(domain(1), domain(2), points);
        x2 = linspace(domain(3), domain(4), points);

        [X1, X2] = ndgrid(x1, x2);
        %
        num = points^2;
        X1 = reshape(X1, num, 1);
        X2 = reshape(X2, num, 1);

        gridX = [X1 X2];

        nmax_needed = num;

        BCs = Create_BCs(ndof,targetPDF,BC_free_sym);

        action = zeros([1, nmax_needed]);
        
        fprintf('%d | tf = %.5f \n', jj,tf)
        parfor_progress(nmax_needed);
        parfor r = 1:nmax_needed
    %         disp(r)

            BVX = gridX(r,:).';

            solinit = bvpinit(tt, zeros(2*ord*ndof, 1));
            % Solve nonlinear BVP
            sol = bvp4c(@M_ode, @M_bc, solinit, options, ODEs, BCs, BVX);

            XX = deval(sol, tt);

            % Calculate action
            Integ = zeros(1,N);
            for it = 1:N
                Integ(it) = Lagrangian([XX(:,it);tt(it)]);
            end

            action(r) = trapz(tt, Integ);

            parfor_progress;
        end
        parfor_progress(0);
            
        PDF = exp(-action);

        elapsed_time = toc;

        save(['files/WPI_EL_joint_' num2str(targetPDF(1)) '_' num2str(targetPDF(2)) ...
                '_tf' num2str(tf) '_fs' num2str(Fs) '.mat'],...
                'points','domain','PDF','targetPDF','par','Fs','elapsed_time');
    end
end