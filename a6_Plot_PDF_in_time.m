clc 
clear all
close all

ndof = 1;
S0 = 0.3;

NR = 1000;
Tot = 10;
Fs_MC = 100.;

points = 31;

targetPDFs = [33];

caxx = 6.5

% aaa = 0.1;
% tfs = [aaa];
% tff = aaa;
% while tff < 1.4
%     tff = round(1.05*tfs(end),4);
%     tfs = [tfs,tff];
% end
% tfs

tfs = [0.100000000000000,0.105000000000000,0.110300000000000,0.115800000000000,...
0.121600000000000,0.127700000000000,0.134100000000000,0.140800000000000,...
0.147800000000000,0.155200000000000,0.163000000000000,0.171200000000000,...
0.179800000000000,0.188800000000000,0.198200000000000,0.208100000000000,...
0.218500000000000,0.229400000000000,0.240900000000000,0.252900000000000,...
0.265500000000000,0.278800000000000,0.292700000000000,0.307300000000000,...
0.322700000000000,0.338800000000000,0.355700000000000,0.373500000000000,...
0.392200000000000,0.411800000000000,0.432400000000000,0.454000000000000,...
0.476700000000000,0.500500000000000,0.525500000000000,0.551800000000000,...
0.579400000000000,0.608400000000000,0.638800000000000,0.670700000000000,...
0.704200000000000,0.739400000000000,0.776400000000000,0.815200000000000,...
0.856000000000000,0.898800000000000,0.943700000000000,0.990900000000000,...
1.040400000000000,1.092400000000000,1.147000000000000,...
1.264600000000000,1.327800000000000,1.394200000000000,1.463900000000000];

Fss_WPI = floor(50./tfs);
Fss_WPI2 = floor(100./tfs);

% tfs_MC = tfs(1):0.01:tfs(end);
% tfs_MC = [tfs,1.21,1.22,1.23,1.24,1.25,1.2646,1.3278,1.3942,1.4639];
tfs_MC = tfs;

XX = [];
YY = [];
ZZ = [];

XX2 = [];
YY2 = [];
ZZ2 = [];

MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];
if exist(MC_file, 'file') == 2
    load(MC_file)
    for k = 1:length(tfs_MC)
        tf = tfs_MC(k);
        n_tf = round(tf*Fs_MC);

        for kk = 1:length(targetPDFs)
            j = targetPDFs(kk);

            X = squeeze(Z(j,n_tf,:));

            [pdf_X_MC,xmesh] = ksdensity(X,'NumPoints',1000);
            XX = [XX;xmesh];
            ZZ = [ZZ;pdf_X_MC];
            YY = [YY;tf*ones(1,length(xmesh))];
        end
    end
end

lgd = {};
for k = 1:length(tfs)
    tf = tfs(k);
    Fs_WPI = Fss_WPI(k);
    Fs_WPI2 = Fss_WPI2(k);
    n_tf = round(tf*Fs_MC);
    
    for kk = 1:length(targetPDFs)
        j = targetPDFs(kk);
        
        WPI_EL_mrg_file = ['WPI_EL/files/WPI_EL_PDF' num2str(j) '_tf' num2str(tf)...
            '_fs' num2str(Fs_WPI) '.mat'];
        WPI_EL_mrg_file2 = ['WPI_EL/files/WPI_EL_PDF' num2str(j) '_tf' num2str(tf)...
            '_fs' num2str(Fs_WPI2) '.mat'];


        % % ---------------------- Path Integral -----------------------------------

        if exist(WPI_EL_mrg_file, 'file') == 2 
            
            load(WPI_EL_mrg_file)

            x = linspace(domain(2*kk-1), domain(2*kk), points);
            xmesh = linspace(domain(2*kk-1), domain(2*kk), 1000);

%             PDF = PDFs{kk};

            total = trapz(x,PDF);
            pdf_x_WPI = PDF/total;
            
            pdf_x_WPIintrp = interp1(x,pdf_x_WPI,xmesh,'pchip');
            
            XX2 = [XX2;xmesh];
            ZZ2 = [ZZ2;pdf_x_WPIintrp];
            YY2 = [YY2;tf*ones(1,length(xmesh))];
            
        elseif  exist(WPI_EL_mrg_file2, 'file') == 2 
            
            load(WPI_EL_mrg_file2)

            x = linspace(domain(2*kk-1), domain(2*kk), points);
            

%             PDF = PDFs{kk};

            total = trapz(x,PDF);
            pdf_x_WPI = PDF/total;
            
            pdf_x_WPIintrp = interp1(x,pdf_x_WPI,xmesh,'pchip');
            
            XX2 = [XX2;xmesh];
            ZZ2 = [ZZ2;pdf_x_WPIintrp];
            YY2 = [YY2;tf*ones(1,length(xmesh))];
            
        end
    end
end

width = 0.3;
height = 0.2;

sz = size(ZZ);
ZZ = [zeros(sz(1),1),ZZ,zeros(sz(1),1)];
YY = [YY(:,1),YY,YY(:,1)];
XX = [ones(sz(1),1),XX,-ones(sz(1),1)];
ZZ2 = [zeros(sz(1),1),ZZ2,zeros(sz(1),1)];
YY2 = [YY2(:,1),YY2,YY2(:,1)];
XX2 = [ones(sz(1),1),XX2,-ones(sz(1),1)];

figure(1);
surf(YY.',XX.',ZZ.','EdgeColor','none','FaceColor','interp')
xlabel('$$t$$','Interpreter','latex','FontSize',40)
ylabel('$$x_{33}$$','Interpreter','latex','FontSize',40)
% xlim([tfs(1),tfs(end)])
% ylim([xmesh(1),xmesh(end)])
set(gca,'fontsize',30)
colorbar()
caxis([0,caxx])
view(2)
set(gca, 'units', 'normalized'); 
Tight = get(gca, 'TightInset');  
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.1 1-Tight(2)-Tight(4)-0.03]; 
set(gca, 'Position', NewPos);
xlim([min(YY(:)),max(YY(:))])
ylim([-0.7,0.7])
set(gca, 'XTick', [0.1,0.4,0.7,1,1.3]);
print('nanoMC','-depsc','-r200')
% print('nanoMC','-dpng','-r1000')

figure(2);
surf(YY2.',XX2.',ZZ2.','EdgeColor','none','FaceColor','interp')
xlabel('$$t$$','Interpreter','latex','FontSize',40)
ylabel('$$x_{33}$$','Interpreter','latex','FontSize',40)
% xlim([tfs(1),tfs(end)])
% ylim([xmesh(1),xmesh(end)])
set(gca,'fontsize',30)
colorbar()
caxis([0,caxx])
view(2)
set(gca, 'units', 'normalized'); 
Tight = get(gca, 'TightInset');  
NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3)-0.1 1-Tight(2)-Tight(4)-0.03];
set(gca, 'Position', NewPos);
xlim([min(YY2(:)),max(YY2(:))])
ylim([-0.7,0.7])
set(gca, 'XTick', [0.1,0.4,0.7,1,1.3]);
print('nanoWPI','-depsc','-r200')
% print('nanoWPI','-dpng','-r1000')


