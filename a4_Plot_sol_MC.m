clc 
clear all
close all

NR = 1000;
Tot = 10.;
Fs_MC = 100.;


MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];

LinT_file = ['MC/files/EqLin_tot' num2str(Tot) '_fs' num2str(Fs_MC)  '.mat'];
  
% figure(1,'Position', [100, 100, 550, 550])
if exist(MC_file, 'file') == 2  
    load(MC_file);
    z1 = squeeze(Z(1,:,:));
    z2 = squeeze(Z(2,:,:));
    z3 = squeeze(Z(3,:,:));
    z4 = squeeze(Z(4,:,:));
    % Time vector (Dimensional)
    [ttMC, dt, nt] = Dim_time( Fs,Tot);

    N = length(z1(1,:,1));
    m_z1 = zeros(N,1);
    m_z2 = zeros(N,1);
    m_z3 = zeros(N,1);
    m_z4 = zeros(N,1);

    var_z1 = zeros(N,1);
    var_z2 = zeros(N,1);
    cov_z1z2 = zeros(N,1);
    var_z3 = zeros(N,1);
    var_z4 = zeros(N,1);
    cov_z3z4 = zeros(N,1);
    for it = 1:N
        m_z1(it) = mean(z1(it,:));
        m_z2(it) = mean(z2(it,:));
        m_z3(it) = mean(z3(it,:));
        m_z4(it) = mean(z4(it,:));

        tmp = cov(z1(it,:),z2(it,:));
        var_z1(it)   = tmp(1,1);
        var_z2(it)   = tmp(2,2);
        cov_z1z2(it) = tmp(1,2);
        tmp_dot = cov(z3(it,:),z4(it,:));
        var_z3(it)   = tmp_dot(1,1);
        var_z4(it)   = tmp_dot(2,2);
        cov_z3z4(it) = tmp_dot(1,2);
    end
    
    figure(1)
    subplot(2,1,1)
    plot(ttMC,m_z1,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
    title('$E(x_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttMC,m_z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$E(x_2)$','interpreter','latex')

    figure(2)
    subplot(2,1,1)
    plot(ttMC,m_z3,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
    title('$E(\dot{x}_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttMC,m_z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$E(\dot{x}_2)$','interpreter','latex')

    figure(3)
    subplot(2,1,1)
    plot(ttMC,var_z1,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
    title('$Var(x_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttMC,var_z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(x_2)$','interpreter','latex')

    figure(4)
    subplot(2,1,1)
    plot(ttMC,var_z3,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(\dot{x}_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttMC,var_z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(\dot{x}_2)$','interpreter','latex')

    figure(5)
    subplot(2,1,1)
    plot(ttMC,cov_z1z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Cov(x_1,x_2)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttMC,cov_z3z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Cov(\dot{x}_1,\dot{x}_2)$','interpreter','latex')
end

if exist(LinT_file, 'file') == 2  
    load(LinT_file)
    ttLin = tt;

    figure(1)
    subplot(2,1,1)
    plot(ttLin,mean_Z(1,:),'color',0.0*[1 1 1],'Linewidth',1);
    title('$E(x_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttLin,mean_Z(2,:),'color',0.0*[1 1 1],'Linewidth',1);
    title('$E(x_2)$','interpreter','latex')
    
    figure(2)
    subplot(2,1,1)
    plot(ttLin,mean_Z(3,:),'color',0.0*[1 1 1],'Linewidth',1);
    title('$E(\dot{x}_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttLin,mean_Z(4,:),'color',0.0*[1 1 1],'Linewidth',1);
    title('$E(\dot{x}_2)$','interpreter','latex')

    figure(3)
    subplot(2,1,1)
    plot(ttLin,squeeze(Cov_Z(1,1,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Var(x_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttLin,squeeze(Cov_Z(2,2,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Var(x_2)$','interpreter','latex')

    figure(4)
    subplot(2,1,1)
    plot(ttLin,squeeze(Cov_Z(3,3,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Var(\dot{x}_1)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttLin,squeeze(Cov_Z(4,4,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Var(\dot{x}_2)$','interpreter','latex')

    figure(5)
    subplot(2,1,1)
    plot(ttLin,squeeze(Cov_Z(1,2,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Cov(x_1,x_2)$','interpreter','latex')
    subplot(2,1,2)
    plot(ttLin,squeeze(Cov_Z(3,4,:)),'color',0.0*[1 1 1],'Linewidth',1);
    title('$Cov(\dot{x}_1,\dot{x}_2)$','interpreter','latex')

end
% squeeze(V(1,1,:))

% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(ttMC,m_z1,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
% plot(ttLin,mean_Z(1,:),'color',0.0*[1 1 1],'Linewidth',1);
% title('$E(x_1)$','interpreter','latex')
% subplot(2,1,2)
% plot(ttMC,m_z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,mean_Z(2,:),'color',0.0*[1 1 1],'Linewidth',1);
% title('$E(x_2)$','interpreter','latex')
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(ttMC,m_z3,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
% plot(ttLin,mean_Z(3,:),'color',0.0*[1 1 1],'Linewidth',1);
% title('$E(\dot{x}_1)$','interpreter','latex')
% subplot(2,1,2)
% plot(ttMC,m_z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,mean_Z(4,:),'color',0.0*[1 1 1],'Linewidth',1);
% title('$E(\dot{x}_2)$','interpreter','latex')
% 
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(ttMC,var_z1,'color',0.5*[1 1 1],'Linewidth',1.5);hold on
% plot(ttLin,squeeze(Cov_Z(1,1,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Var(x_1)$','interpreter','latex')
% subplot(2,1,2)
% plot(ttMC,var_z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,squeeze(Cov_Z(2,2,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Var(x_2)$','interpreter','latex')
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(ttMC,var_z3,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,squeeze(Cov_Z(3,3,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Var(\dot{x}_1)$','interpreter','latex')
% subplot(2,1,2)
% plot(ttMC,var_z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,squeeze(Cov_Z(4,4,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Var(\dot{x}_2)$','interpreter','latex')
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(ttMC,cov_z1z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,squeeze(Cov_Z(1,2,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Cov(x_1,x_2)$','interpreter','latex')
% subplot(2,1,2)
% plot(ttMC,cov_z3z4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% plot(ttLin,squeeze(Cov_Z(3,4,:)),'color',0.0*[1 1 1],'Linewidth',1);
% title('$Cov(\dot{x}_1,\dot{x}_2)$','interpreter','latex')
% 
% 




% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(tt,var_z1,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(1,1)) real(Cov(1,1))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Var(y_1)$','interpreter','latex')
% legend('MC500 non linear','Stationary - Freq. Domain','Location','southeast')
% subplot(2,1,2)
% plot(tt,var_z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(2,2)) real(Cov(2,2))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Var(y_2)$','interpreter','latex')
% legend('MC500 non linear','Stationary - Freq. Domain','Location','southeast')
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(tt,var_Z1,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(3,3)) real(Cov(3,3))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Var(\dot{y}_1)$','interpreter','latex')
% legend('MC500 non linear','Stationary - Freq. Domain','Location','southeast')
% subplot(2,1,2)
% plot(tt,var_Z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(4,4)) real(Cov(4,4))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Var(\dot{y}_2)$','interpreter','latex')
% legend('MC500 non linear','Stationary - Freq. Domain','Location','southeast')
% 
% figure('Position', [100, 100, 550, 550])
% subplot(2,1,1)
% plot(tt,cov_z1z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(1,2)) real(Cov(1,2))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Cov(y_1,y_2)$','interpreter','latex')
% legend('MC500 non linear' ,'Stationary - Freq. Domain','Location','southeast')
% subplot(2,1,2)
% plot(tt,cov_Z1Z2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
% % plot([1 tt(end)],[real(Cov(3,4)) real(Cov(3,4))],'--','color',0.0*[1 1 1],'Linewidth',1.5)
% xlabel('time (sec)')
% title('$Cov(\dot{y}_1,\dot{y}_2)$','interpreter','latex')
% legend('MC500 non linear','Stationary - Freq. Domain','Location','southeast')