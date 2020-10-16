% load('MC500_CP_non_lin');
function [] = plot_mean_var( z1,z2,z3,z4,dt )
    N = length(z1(:,1));
    tt = (0:N-1)*dt;       % Time vector
    
    m_1 = zeros(N,1);
    m_2 = zeros(N,1);
    m_3 = zeros(N,1);
    m_4 = zeros(N,1);
    
    var_1 = zeros(N,1);
    var_2 = zeros(N,1);
    var_3 = zeros(N,1);
    var_4 = zeros(N,1);
    
    for it = 1:N
        m_1(it) = mean(z1(it,:));
        m_2(it) = mean(z2(it,:));
        m_3(it) = mean(z3(it,:));
        m_4(it) = mean(z4(it,:));
       
        var_1(it)     = cov(z1(it,:));
        var_2(it)     = cov(z2(it,:));
        var_3(it)     = cov(z3(it,:));
        var_4(it)     = cov(z4(it,:));
        
    end
    
    figure('Position', [100, 100, 550, 550])
    subplot(4,1,1)
    plot(tt,m_1,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
%     plot(tt,ones(N)*Lm_x,'color',0.0*[1 1 1],'Linewidth',1.5); hold on
    title('$E(x_1)$','interpreter','latex')
    subplot(4,1,2)
    plot(tt,m_2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$E(x_2)$','interpreter','latex') 
    subplot(4,1,3)
    plot(tt,m_3,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$E(\dot{x_1})$','interpreter','latex')
    subplot(4,1,4)
    plot(tt,m_4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$E(\dot{x_2})$','interpreter','latex')
    xlabel('time (sec)')
    
    
    figure('Position', [100, 100, 550, 550])
    subplot(4,1,1)
    plot(tt,var_1,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
%     plot(tt,ones(N)*Lvar_x,'color',0.0*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(x_1)$','interpreter','latex')
    subplot(4,1,2)
    plot(tt,var_2,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(x_1)$','interpreter','latex')
    subplot(4,1,3)
    plot(tt,var_3,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(\dot{x_1})$','interpreter','latex')
    subplot(4,1,4)
    plot(tt,var_4,'color',0.5*[1 1 1],'Linewidth',1.5); hold on
    title('$Var(\dot{x_2})$','interpreter','latex')
    xlabel('time (sec)')
end

