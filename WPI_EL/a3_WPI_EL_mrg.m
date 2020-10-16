clc
clear all
close all
tic

mydir  = pwd;
idcs   = strfind(mydir,'\');
newdir = mydir(1:idcs(end)-1);
load([newdir '/par.mat'])

MC_file = [newdir '/MC/files/MC45_' num2str(1000) '_tot' num2str(10) '_fs' num2str(100) '.mat'];

ndof = par.ndof;
ord  = 2;
ti = 0;

% tfs = [3.];
% tfs = [0.5 1.0 1.5 2.0 2.5 3.0 3.5 4.0 4.5 5.0];

% aaa = 1.3942;
% tfs = [aaa];
% tff = aaa;
% while tff < 1.7
%     tff = round(1.05*tfs(end),4);
%     tfs = [tfs,tff];
% end
tfs = [1.21,1.22,1.23,1.24,1.25];
disp(size(tfs))

% Time parameters
Fss = floor(100./tfs);

points = 31 ;
targetPDFs = [33];

b3_functional_minimization_EL_marginal( MC_file, points, ord, targetPDFs, par, ti, tfs, Fss);

% beep