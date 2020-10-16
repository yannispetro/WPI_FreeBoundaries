clc 
clear all
close all

ndof = 100;

NR = 1000;
Tot = 10.;
Fs_MC = 100.;

tf = 0.5;
n_tf = round(tf*Fs_MC);

Fs_WPI = 100.;
points = 31;

targetPDFs = [1];

MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];

WPI_EL_mrg_file = ['WPI_EL/files/WPI_EL_tf' num2str(tf) '_fs' num2str(Fs_WPI) '.mat'];

% --------------------- MC ---------------------------------------
if exist(MC_file, 'file') == 2
    load(MC_file)
end

for kk = 1:length(targetPDFs)
    j = targetPDFs(kk);
    
    if exist(MC_file, 'file') == 2

        X = squeeze(Z(j,n_tf,:));

        [pdf_X_MC,xmesh] = ksdensity(X,'NumPoints',100);
        
        figure(kk); plot(xmesh,pdf_X_MC,'LineStyle','--','Linewidth',2); hold on

    end


    % % ---------------------- Path Integral -----------------------------------

    if exist(WPI_EL_mrg_file, 'file') == 2 
        load(WPI_EL_mrg_file)

        x = linspace(domain(2*kk-1), domain(2*kk), points);
        
        PDF = PDFs{kk};

        total = trapz(x,PDF);
        pdf_x_WPI = PDF/total;

        figure(kk); plot(x,pdf_x_WPI, '-x','LineStyle','-','Linewidth',1,'MarkerSize',10); hold on
%         xlim([x(1),x(end)])

        figure(kk);
        if j <= ndof
            var = ['{x}_{' num2str(j) '}'];
        else
            var = ['\dot{{x}}_{' num2str(j-ndof) '}'];
        end       
        xlabel(['$$' var '$$'],'Interpreter','latex','FontSize',30)
        title(['PDF of $$' var '$$'],'Interpreter','latex','FontSize',20)
        legend(['t = ' num2str(tf) ' s - MC' num2str(NR)],['t = ' num2str(tf) ' s - WPI'])
    end
end
