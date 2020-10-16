function domain = c1_Get_bounds(MC_file,targetPDFs,tf,tol)
    load(MC_file);
    
%     s = size(Z);
%     nsdof = s(1);
    n_dom = length(targetPDFs);
    
    n_tf = round(tf*Fs);
    
    [ tt,dt,N ] = Dim_time( Fs,Tot );
    
    X = squeeze(Z(:,n_tf,:));
    boundzz = [];
    for i = 1:n_dom
        j = targetPDFs(i);
        x = X(j,:);
        [pdf_x_MC,xmesh] = ksdensity(x,'NumPoints',1000);
        m_pdf = max(pdf_x_MC);
        x_n = xmesh( pdf_x_MC(:) > tol*m_pdf );
        xl = x_n(1);
        xu = x_n(end);
        boundzz = [boundzz,xl,xu];  
%         xl = round(10*x_n(1))./10;
%         if xl == 0
%             xl = x_n(1);
%         end
%         xu = round(10*x_n(end))./10;
%         if xu == 0
%             xu = x_n(end);
%         end
%         xlu = max(abs(xl),abs(xu));
%         boundzz = [boundzz,xl,xu];  
    end
    domain = boundzz;
    
    
%     its = 5:10:n_tf;
%     if its(end) < n_tf
%         its2 = [its,n_tf];
%     else
%         its2 = its;
%     end
% 
%     dom_t = zeros(n_dom,length(its2));
%     for jj = 1:length(its2)
%         it = its2(jj);
%         n_t = round(tt(it)*Fs);
%         X = squeeze(Z(:,n_t,:));
%         boundzz = [];
%         for i = 1:n_dom
%             j = targetPDFs(i);
%             x = X(j,:);
%             [pdf_x_MC,xmesh] = ksdensity(x,'NumPoints',1000);
%             x_n = xmesh( pdf_x_MC(:) > 1e-4 ) ;
%             xl = round(10*x_n(1))./10;
%             xu = round(10*x_n(end))./10;
%             xlu = max(abs(xl),abs(xu));
%             boundzz = [boundzz,xlu];  
%         end
%         dom_t(:,jj) = boundzz;
%     end
%     domain = boundzz;
%     
%     tp = [0,tt(its2)];
%     S = [zeros(n_dom,1),dom_t(:,1:end)];
%     plot( tp,S );
% 
%     drawnow;
end

