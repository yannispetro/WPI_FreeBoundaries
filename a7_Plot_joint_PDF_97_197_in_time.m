clc 
clear all
close all

ndof = 10;

NR = 1000;
Tot = 10.;
Fs_MC = 100.;

points = 21;
pointsMC = 51;
nkde2 = 2^8;

targetPDF = [97,197];

% ts1 = 0.025*(1:2:23);
% ts2 = 0.05*(1:1:19);
% % tfs = [0.015,0.02,ts1,ts2];
% % Fss = [200,200,200*ones(1,length(ts1)),100*ones(1,length(ts2))];
% % [tfs, tf_order] = sort(tfs);
% % Fss = Fss(tf_order);
% 
% tfs = [0.025,ts2];
% Fss = [200,100*ones(1,length(ts2))];

tfs = (0.025:0.05:1.075);
% tfs = (0.05:0.05:0.9);
Fss = 200*ones(1,length(tfs));


plevels = [0.01,0.32,0.6];
alphas  = [0.3,0.3,1];
side_alpha = 0.9;
cax = [-0.02 0.07];

j1 = targetPDF(1);
j2 = targetPDF(2);
% --------------------- MC ---------------------------------------
MC_file = ['MC/files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs_MC) '.mat'];
load(MC_file)
[tt, dt, NtMC] = Dim_time( Fs,Tot );

XV = zeros(length(tfs),points,points);
YV = zeros(length(tfs),points,points);
ZV = zeros(length(tfs),points,points);
V = zeros(length(tfs),points,points);
XVmc = zeros(length(tfs),pointsMC,pointsMC);
YVmc = zeros(length(tfs),pointsMC,pointsMC);
ZVmc = zeros(length(tfs),pointsMC,pointsMC);
Vmc = zeros(length(tfs),pointsMC,pointsMC);
for k = 1:length(tfs)
    Fs_WPI = Fss(k);
    tf = tfs(k);
    
    WPI_EL_joint_file = ['WPI_EL/files/WPI_EL_joint_' num2str(j1) '_' num2str(j2) ...
                         '_tf' num2str(tf) '_fs' num2str(Fs_WPI) '.mat'];

    load(WPI_EL_joint_file) 
    
    n_tf = round(tf*Fs_MC);
    x1 = squeeze(Z(j1,n_tf,:));
    x2 = squeeze(Z(j2,n_tf,:));

    [~,density,Xmesh,Ymesh]=kde2d([x1,x2],nkde2,[domain(1), domain(3)], [domain(2), domain(4)]);
    XX1mc = linspace(domain(1), domain(2), pointsMC);
    XX2mc = linspace(domain(3), domain(4), pointsMC);
    for i = 1:pointsMC
        XVmc(k,i,:) = XX1mc;
        YVmc(k,:,i) = XX2mc;
    end
    ZVmc(k,:,:) = tf;
    
    [Xinterp,Yinterp] = meshgrid(XX1mc,XX2mc);
    PDF2Dmc = interp2(Xmesh,Ymesh,density,Xinterp,Yinterp);
    
    I2D = trapz(XX2mc,trapz(XX1mc,PDF2Dmc,2));
    Vmc(k,:,:) = PDF2Dmc/I2D;

    
    
% ------------ Path Integral -------------------------
    XX1 = linspace(domain(1), domain(2), points);
    XX2 = linspace(domain(3), domain(4), points);
    for i = 1:points
        XV(k,i,:) = XX1;
        YV(k,:,i) = XX2;
    end
    ZV(k,:,:) = tf;
    
    PDF2D = reshape(PDF, points,points);
    PDF2Dcor = PDF2D;
    for i = 2:length(PDF2D)-1
        for j = 2:length(PDF2D)-1
            val = ( PDF2Dcor(i-1,j) + PDF2Dcor(i,j-1) + ...
                    PDF2Dcor(i+1,j) + PDF2Dcor(i,j+1) + ...
                    PDF2Dcor(i-1,j-1) + PDF2Dcor(i+1,j-1) + ...
                    PDF2Dcor(i+1,j+1) + PDF2Dcor(i-1,j+1) )/8;
            if PDF2Dcor(i,j)/val < 0.8
                PDF2Dcor(i,j) = val;
            end
        end
    end
    I2D = trapz(XX2,trapz(XX1,PDF2Dcor,2));
    V(k,:,:) = PDF2Dcor/I2D;
    
end

figMC = figure('Position',  [440   378   560   390]);
axMC = gca;
set(axMC,'Ydir','reverse')

iso1mc = patch(isosurface(ZVmc,XVmc,YVmc,Vmc,plevels(1)));
iso2mc = patch(isosurface(ZVmc,XVmc,YVmc,Vmc,plevels(2)));
iso3mc = patch(isosurface(ZVmc,XVmc,YVmc,Vmc,plevels(3)));

iso1mc.EdgeColor = 'none';
iso2mc.EdgeColor = 'none';
iso3mc.EdgeColor = 'none';

iso1mc.FaceAlpha = alphas(1);
iso2mc.FaceAlpha = alphas(2);
iso3mc.FaceAlpha = alphas(3);


Nmc = 10*pointsMC;
XXXmc = linspace(axMC.YLim(1),axMC.YLim(2),Nmc);
YYYmc = linspace(axMC.ZLim(1),axMC.ZLim(2),Nmc);
V_Xtmc = zeros(Nmc,length(tfs),Nmc);
V_Ytmc = zeros(Nmc,length(tfs),Nmc);
for k = 1:length(tfs)
    PDF2D = squeeze(Vmc(k,:,:));
    X_orig = squeeze(XVmc(k,1,:)).';
    Y_orig = squeeze(YVmc(k,:,1));

    PDF_Xt_orig = Nmc^(-1)*sum(PDF2D,1);
    PDF_Yt_orig = Nmc^(-1)*sum(PDF2D,2).';
%     PDF_Xt_orig = PDF_Xt_orig/trapz(X_orig,PDF_Xt_orig);
%     PDF_Yt_orig = PDF_Yt_orig/trapz(Y_orig,PDF_Yt_orig);
%     
    PDF_Xt = interp1([XXXmc(1),X_orig,XXXmc(end)],[0,PDF_Xt_orig,0],XXXmc);
    PDF_Yt = interp1([YYYmc(1),Y_orig,YYYmc(end)],[0,PDF_Yt_orig,0],YYYmc);
    
    for i = 1:Nmc
        V_Xtmc(:,k,i) = PDF_Xt;
        V_Ytmc(i,k,:) = PDF_Yt;
    end
end

[X,Y,Z] = meshgrid(tfs,XXXmc,YYYmc);
hold on
hxmc = slice(X,Y,Z,V_Xtmc,[],[],axMC.ZLim(1));
set(hxmc,'edgecolor','none')
set(hxmc,'facecolor','interp')
set(hxmc,'facealpha',side_alpha)

hold on
hymc = slice(X,Y,Z,V_Ytmc,[],axMC.YLim(1),[]);
set(hymc,'edgecolor','none')
set(hymc,'facecolor','interp')
set(hymc,'facealpha',side_alpha)

caxis(cax)

pbaspect([3 1 1])
%set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',15)
view(3); 
axis tight
camlight 
lighting gouraud

xlh = xlabel('$$t$$','Interpreter','latex','FontSize',20);
ylh = ylabel('$$x_{97}$$','Interpreter','latex','FontSize',20);
zlh = zlabel('$\dot{x}_{97}$','Interpreter','latex','FontSize',20);

xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.6);
xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.4);

ylh.Position(1) = ylh.Position(1) + abs(ylh.Position(1) * 1.5);
ylh.Position(2) = ylh.Position(2) - abs(ylh.Position(2) * 8);

zlh.Position(1) = zlh.Position(1) + abs(zlh.Position(1) * 1);
zlh.Position(2) = zlh.Position(2) - abs(zlh.Position(2) * 0.01);

outerpos = axMC.OuterPosition;
ti = axMC.TightInset; 
left = outerpos(1) + ti(1)*1.7;
bottom = outerpos(2) - ti(2)*9.3;
ax_width = outerpos(3) - ti(1)*2 - ti(3);
ax_height = outerpos(4) + ti(2)*20 - ti(4);
axMC.Position = [left bottom ax_width ax_height];

box on


% ------------ Path Integral -------------------------
figWPI = figure('Position',  [440   378   560   390]);
ax = gca;
set(ax,'Ydir','reverse')

iso1 = patch(isosurface(ZV,XV,YV,V,plevels(1)));
iso2 = patch(isosurface(ZV,XV,YV,V,plevels(2)));
iso3 = patch(isosurface(ZV,XV,YV,V,plevels(3)));

iso1.EdgeColor = 'none';
iso2.EdgeColor = 'none';
iso3.EdgeColor = 'none';

iso1.FaceAlpha = alphas(1);
iso2.FaceAlpha = alphas(2);
iso3.FaceAlpha = alphas(3);


N = 10*points;
ax.YLim = ax.YLim*1.01;
ax.ZLim = ax.ZLim*1.01;
XXX = linspace(ax.YLim(1),ax.YLim(2),N);
YYY = linspace(ax.ZLim(1),ax.ZLim(2),N);
V_Xt = zeros(N,length(tfs),N);
V_Yt = zeros(N,length(tfs),N);
for k = 1:length(tfs)
    PDF2D = squeeze(V(k,:,:));
    X_orig = squeeze(XV(k,1,:)).';
    Y_orig = squeeze(YV(k,:,1));

    PDF_Xt_orig = N^(-1)*sum(PDF2D,1);
    PDF_Yt_orig = N^(-1)*sum(PDF2D,2).';
%     PDF_Xt_orig = PDF_Xt_orig/trapz(X_orig,PDF_Xt_orig);
%     PDF_Yt_orig = PDF_Yt_orig/trapz(Y_orig,PDF_Yt_orig);
%     
    PDF_Xt = interp1([XXX(1),X_orig,XXX(end)],[0,PDF_Xt_orig,0],XXX);
    PDF_Yt = interp1([YYY(1),Y_orig,YYY(end)],[0,PDF_Yt_orig,0],YYY);
    for i = 1:N
        V_Xt(:,k,i) = PDF_Xt;
        V_Yt(i,k,:) = PDF_Yt;
    end
end

[X,Y,Z] = meshgrid(tfs,XXX,YYY);
hold on
hx = slice(X,Y,Z,V_Xt,[],[],ax.ZLim(1));
set(hx,'edgecolor','none')
set(hx,'facecolor','interp')
set(hx,'facealpha',side_alpha)

hold on
hy = slice(X,Y,Z,V_Yt,[],ax.YLim(1),[]);
set(hy,'edgecolor','none')
set(hy,'facecolor','interp')
set(hy,'facealpha',side_alpha)

caxis(cax)

pbaspect([3 1 1])
% set(gcf,'renderer','zbuffer') 
set(gca,'fontsize',15)
view(3); 
axis tight
camlight 
lighting gouraud

xlh = xlabel('$$t$$','Interpreter','latex','FontSize',20);
ylh = ylabel('$$x_{97}$$','Interpreter','latex','FontSize',20);
zlh = zlabel('$\dot{x}_{97}$','Interpreter','latex','FontSize',20);


xlh.Position(1) = xlh.Position(1) + abs(xlh.Position(1) * 0.6);
xlh.Position(2) = xlh.Position(2) - abs(xlh.Position(2) * 0.4);

ylh.Position(1) = ylh.Position(1) + abs(ylh.Position(1) * 1.2);
ylh.Position(2) = ylh.Position(2) - abs(ylh.Position(2) * 15);

zlh.Position(1) = zlh.Position(1) + abs(zlh.Position(1) * 1);
zlh.Position(2) = zlh.Position(2) - abs(zlh.Position(2) * 0.01);

outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1)*1.7;
bottom = outerpos(2) - ti(2)*9.3;
ax_width = outerpos(3) - ti(1)*2 - ti(3);
ax_height = outerpos(4) + ti(2)*20 - ti(4);
ax.Position = [left bottom ax_width ax_height];

box on

% ----------- Both MC and WPI -------------------
clm = colormap;

iso1mc.FaceColor = clm(25,:);
iso2mc.FaceColor = clm(53,:);
iso3mc.FaceColor = clm(end,:);

iso1.FaceColor = clm(25,:);
iso2.FaceColor = clm(53,:);
iso3.FaceColor = clm(end,:);

print(figMC,'nanoJointMCS_time','-depsc','-r200')
print(figWPI,'nanoJointWPI_time','-depsc','-r200')

% print(figMC,'nanoJointMCS_time','-dpng','-r1000')
% print(figWPI,'nanoJointWPI_time','-dpng','-r1000')
