%% Matthew Ning
% Adapted from "main.m" adapted by Jasmine
% This calculates condition-averaged PLI ONLY

close all
clear 
clc
tic

global FN BN FS1 BS1 FS2 BS2;
FN = 1;
BN = 2;
FS1 = 3;
BS1 = 4;
FS2 = 5;
BS2 = 6;
% Load raw data
subjINFO = 'ASA067_XX_8-14_DATA.mat';

home = '~/Documents/RotationSpring2018_2';
load(subjINFO)

CH = size(ERP.erp,1);
ts = size(ERP.erp,2);
trials = size(ERP.erp,3);
phit = zeros(CH,CH,trials);
PLIt =phit;

[PLI, pliTrial] = calcSync_Stam(ERP.erp,CH,SCORE.F_N,SCORE.B_N,SCORE.F_S1,SCORE.B_S1,SCORE.F_S2,SCORE.B_S2);
[freq, PLV, plvTrial, PLIV, pliTrialV, WPLI,wpliTrial, Coherency] = ...
    calcSync_Vinck(ERP.erp,SCORE.F_N,SCORE.B_N,SCORE.F_S1,SCORE.B_S1,SCORE.F_S2,SCORE.B_S2);

names = {'Fp1'; 'AF7'; 'AF3'; 'F1'; 'F3'; 'F5'; 'F7'; 'FT7'; 'FC5'; 'FC3'; 'FC1';...
    'C1'; 'C3'; 'C5'; 'T7'; 'TP7'; 'CP5'; 'CP3'; 'CP1'; ...
    'P1'; 'P3'; 'P5'; 'P7'; 'P9'; 'P07'; 'P03'; '01'; 'Iz'; '0z';...
    'P0z'; 'Pz';'CPz'; 'Fpz'; 'Fp2'; 'AF8'; 'AF4'; 'Afz'; ...
    'Fz'; 'F2'; 'F4'; 'F6'; 'F8'; 'FT8'; 'FC6'; 'FC4'; 'FC2'; 'FCz'; ...
    'Cz'; 'C2'; 'C4'; 'C6'; 'T8'; 'TP8'; 'CP6'; 'CP4'; 'CP2'; ...
    'P2'; 'P4'; 'P6'; 'P8'; 'P10'; 'P08'; 'P04'; '02'; '1'; '2'; '3';};
save([home,'/data/', strcat(subjINFO(1:9),'_NETWORK_PLI_All.mat')])
%% Calculate Instantaneous Phase (if it returns wrapped (-pi to pi), unwrap it (-inf to inf))
% Leave it wrapped, radian

% In all metrics, when computing expected value, assume relative phase is
% uniformly distributed on circular distribution. So E[X(phi)] = mean(phi)

%%  Plot for proof of concept
chi = 12;
chj = 24;
trial = 4; 
t = linspace(-2.5,4.5,1792);

figure(1);
subplot(2,1,1)
plot(t,squeeze(ERP.erp(chi,:,trial)));
hold on;
plot(t,squeeze(ERP.erp(chj,:,trial)));
xlabel('Time [ms]');
ylabel('Amplitude [mV]');
legend(['Channel: ',num2str(chi)],['Channel: ',num2str(chj)]);
title(['EEG (PLV: ',num2str(pliTrial(chi,chj,trial)),')']);

% subplot(2,1,2)
% hold on
% plot(t,squeeze(mean(dPHI_FN(chi,chj,:),3)))
% plot(t,squeeze(mean(dPHI_BN(chi,chj,:),3)))

%% Contrast plots (No supertarget) 
figure(2);
imagesc(PLI(:,:,FN));
colormap;
title(strcat('Stam Subj ',subjINFO(5:6),': PLI Focus (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(3);
imagesc(PLI(:,:,BN));
colormap;
title(strcat('Stam Subj ',subjINFO(5:6),': PLI Broad (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(4);
imagesc(PLIV(:,:,FN));
colormap;
title(strcat('Vinck Subj ',subjINFO(5:6),': PLI Focus (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(5);
imagesc(PLIV(:,:,BN));
colormap;
title(strcat('Vinck Subj ',subjINFO(5:6),': PLI Broad (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(6);
imagesc(PLV(:,:,FN));
colormap;
title(strcat('Subj ',subjINFO(5:6),': PLV Focus (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(7);
imagesc(PLV(:,:,BN));
colormap;
title(strcat('Subj ',subjINFO(5:6),': PLV Broad (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(8);
imagesc(WPLI(:,:,FN));
colormap;
title(strcat('Subj ',subjINFO(5:6),': WPLI Focus (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

figure(9);
imagesc(WPLI(:,:,BN));
colormap;
title(strcat('Subj ',subjINFO(5:6),': WPLI Broad (No Supertarget)'));
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);
colorbar

toc