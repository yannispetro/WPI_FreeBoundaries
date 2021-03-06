clc 
clear all
close all

ndof = 10;

NR = 1000;
Tot = 10.;
Fs_MC = 100.;

tf = 0.625;
n_tf = round(tf*Fs_MC);

Fs_WPI = 200.;
points = 21;
targetPDF = [97,197];
% targetPDF = [96,97];
nkde2 = 2^8;

j1 = targetPDF(1);
j2 = targetPDF(2);
MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];

WPI_EL_joint_file = ['WPI_EL/files/WPI_EL_joint_' num2str(j1) '_' num2str(j2) ...
                     '_tf' num2str(tf) '_fs' num2str(Fs_WPI) '.mat'];

load(WPI_EL_joint_file)
load(MC_file)
% --------------------- MC ---------------------------------------
[tt, dt, N] = Dim_time( Fs,Tot );
x1 = squeeze(Z(j1,n_tf,:));
x2 = squeeze(Z(j2,n_tf,:));

[~,density,Xmesh,Ymesh]=kde2d([x1,x2],nkde2,[domain(1), domain(3)], [domain(2), domain(4)]);

cax = max(density(:));

figure();
surf(Xmesh,Ymesh,density,'FaceColor','interp','LineStyle','none')
xlabel('$$x_{100}$$','Interpreter','latex','FontSize',40)
ylabel('$$\dot{x}_{100}$$','Interpreter','latex','FontSize',40)
zlabel('$PDF$','Interpreter','latex')
xlim([domain(1),domain(2)]);
ylim([domain(3),domain(4)]);
% xlim([-0.46,0.32]);
% ylim([-8.2,8.0]);
set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',30)
view(2)
colorbar()
caxis([0,cax])
print('nanoJointMCS','-depsc','-r1000')


% % ---------------------- Path Integral -----------------------------------
load(WPI_EL_joint_file);
XX1 = linspace(domain(1), domain(2), points);
XX2 = linspace(domain(3), domain(4), points);
PDF2D = reshape(PDF, points,points);

I2D = trapz(XX2,trapz(XX1,PDF2D,2));
PDF2D = PDF2D/I2D;

figure();
surf(XX1,XX2,PDF2D.','FaceColor','interp','LineStyle','none')
xlabel('$$x_{100}$$','Interpreter','latex','FontSize',40)
ylabel('$$\dot{x}_{100}$$','Interpreter','latex','FontSize',40)
zlabel('$PDF$','Interpreter','latex')
xlim([domain(1),domain(2)]);
ylim([domain(3),domain(4)]);
% xlim([-0.46,0.32]);
% ylim([-8.2,8.0]);
set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',30)
view(2)
colorbar()
caxis([0,cax])
print('nanoJointWPI','-depsc','-r1000')

PDF2Dcor = PDF2D;
np = 0;
for i = 2:length(PDF2D)-1
    for j = 2:length(PDF2D)-1
        val = ( PDF2Dcor(i-1,j) + PDF2Dcor(i,j-1) + ...
                PDF2Dcor(i+1,j) + PDF2Dcor(i,j+1) + ...
                PDF2Dcor(i-1,j-1) + PDF2Dcor(i+1,j-1) + ...
                PDF2Dcor(i+1,j+1) + PDF2Dcor(i-1,j+1) )/8;
        if PDF2Dcor(i,j)/val < 0.8
            np = np + 1;
            PDF2Dcor(i,j) = val;
        end
    end
end
disp(np)

figure();
surf(XX1,XX2,PDF2Dcor.','FaceColor','interp','LineStyle','none')
xlabel('$$x_{100}$$','Interpreter','latex','FontSize',40)
ylabel('$$\dot{x}_{100}$$','Interpreter','latex','FontSize',40)
zlabel('$PDF$','Interpreter','latex')
xlim([domain(1),domain(2)]);
ylim([domain(3),domain(4)]);
% xlim([-0.46,0.32]);
% ylim([-8.2,8.0]);
set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',30)
view(2)
colorbar()
caxis([0,cax])
%print('nanoJointWPI','-depsc','-r1000')
