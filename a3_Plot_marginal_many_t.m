clc 
clear all
close all

ndof = 100;

NR = 1000;
Tot = 10.;
Fs_MC = 100.;

points = 31;

targetPDFs = [33];

% tfs = [0.2 0.9 1.2 1.4];
tfs = [0.7042,1.21,1.22,1.23,1.24];
Fss_WPI = floor(100./tfs);


MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];
if exist(MC_file, 'file') == 2
    load(MC_file)
end

lgd = {};
for k = 1:length(tfs)
    tf = tfs(k);
    Fs_WPI = Fss_WPI(k);
    n_tf = round(tf*Fs_MC);
    
    
    
    for kk = 1:length(targetPDFs)
        j = targetPDFs(kk);
        
        WPI_EL_mrg_file = ['WPI_EL/files/WPI_EL_PDF' num2str(j) '_tf' num2str(tf) '_fs' num2str(Fs_WPI) '.mat'];

        if exist(MC_file, 'file') == 2

            X = squeeze(Z(j,n_tf,:));
            mx = mean(X);
            
%             disp([j,tf,mx])

            [pdf_X_MC,xmesh] = ksdensity(X,'NumPoints',1000);

            figure(kk); plot(xmesh,pdf_X_MC,'LineStyle','--','Linewidth',2); hold on
            lgd{end+1} = ['t = ' num2str(tfs(k)) ' s - MC' num2str(NR)];
            
%             yl = ylim;
%             plot([mx,mx],[0,yl(2)])

        end

        % % ---------------------- Path Integral -----------------------------------

        if exist(WPI_EL_mrg_file, 'file') == 2 
            load(WPI_EL_mrg_file)

            x = linspace(domain(2*kk-1), domain(2*kk), points);

%             PDF = PDFs{kk};

            total = trapz(x,PDF);
            pdf_x_WPI = PDF/total;

            figure(kk); plot(x,pdf_x_WPI, '-x','LineStyle','-','Linewidth',1,'MarkerSize',10); hold on
            lgd{end+1} = ['t = ' num2str(tfs(k)) ' s - WPI'];
%             xlim([x(1),x(end)])

            figure(kk);
            if j <= ndof
                var = ['{x}_{' num2str(j) '}'];
            else
                var = ['\dot{{x}}_{' num2str(j-ndof) '}'];
            end      
            xlabel(['$$' var '$$'],'Interpreter','latex','FontSize',30)
            title(['PDF of $$' var '$$'],'Interpreter','latex','FontSize',20)
            
        end
        legend(lgd)
    end
    
end

