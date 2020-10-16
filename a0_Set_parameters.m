par.ndof = 100;

par.m0 = 1.;
par.c0 = 1.5;
par.k0 = 120;
par.w0 = 2*pi*35.1;

par.d0 = 0.1;

par.e1 = 0.1;
par.e2 = 0.1;

par.Acos = 20;
par.wcos = 3*pi;

par.S0 = 5.;

par.D0 = eye(par.ndof);

save('par.mat','par')

% ORIGINAL VALUES
% par.m0 = 1.;
% par.c0 = 0.002;
% par.k0 = 0.01;
% 
% par.e1 = 1.;
% par.e2 = 3.;
% 
% par.Acos = 10.;
% par.wcos = pi;
% 
% par.S0 = 1.0;