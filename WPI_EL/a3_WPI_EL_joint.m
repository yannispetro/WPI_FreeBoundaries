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

tfs = 0.05*[19:30];

% Time parameters
Fss = 200*ones(length(tfs),1);
% Fss = floor(10./tfs);

points = 21;
targetPDF = [97,197];

tic
b3_functional_minimization_EL_joint( ... 
           MC_file, points, ord, targetPDF, ... 
           par, ti, tfs, Fss);

toc

beep