clc
clear all
close all
tic

mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);

load([newdir '\par.mat'])

ndof = par.ndof;

m0   = par.m0;
c0   = par.c0;
k0   = par.k0;
w0   = par.w0;

d0   = par.d0;

e1   = par.e1;
e2   = par.e2;

Acos = par.Acos;
wcos = par.wcos;

S0   = par.S0;

D0    = par.D0;

M = m0*eye(ndof);
C = c0*eye(ndof);
K = w0*eye(ndof) + k0*full(gallery('tridiag',ndof,-1,2,-1));

G1 = e1*eye(ndof);
G2 = e2*eye(ndof);

% Time vector (Dimensional)
Fs   = 100;  % Sampling frequency      
Tot  = 10;  % Total time  (sec)
[tt, dt, nt] = Dim_time( Fs,Tot);   

% ----- Modulating function ----
At = ones([ndof,nt]); % sin(2.5*pi*tt/tt(end)).^2 + 0.2; % Time modulation

% ----- Deterministic Force ----
Fdet1 = Acos*cos(wcos*tt); % Time modulation
Fdet = [];
for i = 1:ndof
    Fdet = [Fdet;Fdet1];
end

% Number of realizations
NR = 1000;
% Samples = zeros([ndof,nt,NR]);
tspan = [tt(1) tt(end)];
x0 = zeros(1,2*ndof);
Z = zeros([2*ndof,nt,NR]);
parfor_progress(NR);
parfor jj = 1:NR
%     jj
    fsr = zeros(ndof,nt);
    for i = 1:ndof
        fsr(i,:) = wgn(nt,1,db(sqrt(2*pi*S0/(dt))));
    end
    f = (-1)*At.*fsr;
%     F = D0*f;
    F = Fdet + D0*f;

    My = inv(M);
    [t1,response_temp] = ode45(@(t,x) ODE_MC(t,x,My,C,K,d0,G1,G2,tt,F),tspan,x0);
    
    ZZ = interp1(t1,response_temp,tt,'linear');
    Z(:,:,jj) = ZZ.';

    parfor_progress;
end
parfor_progress(0);

toc
elapsed_time = toc;

save(['files/MC45_' num2str(NR) '_tot' num2str(Tot) '_fs' num2str(Fs) '.mat'],...
      'Z','par','Fs','Tot','elapsed_time', '-v7.3');
  
% plot_mean_var( z1,z2,z3,z4,dt )

beep