%% Matthew Ning
% Auditory Neuroscience Lab
% Professor Barbara Shinn-Cunningham
% Supervised by Jasmine Kwasa

close all
clear all
clc

% Load raw data
load('ASA028_XX_8-14_DATA.mat');
CH = size(ERP.erp,1);
ts = size(ERP.erp,2);
trials = size(ERP.erp,3);
phit = zeros(CH,CH,trials);
PLIt =phit;
WPLIt = phit;
Ct = phit;

FS1s = size(ERP.erp(:,:,SCORE.F_S1),3);
FS2s = size(ERP.erp(:,:,SCORE.F_S2),3);
FNs = size(ERP.erp(:,:,SCORE.F_N),3);
BS1s = size(ERP.erp(:,:,SCORE.B_S1),3);
BS2s = size(ERP.erp(:,:,SCORE.B_S2),3);
BNs = size(ERP.erp(:,:,SCORE.B_N),3);

PLVFS1 = zeros(CH,ts,1);
PLVFS2 = PLVFS1;
PLVFN = PLVFS1;
PLVBS1 = PLVFS1;
PLVBS2 = PLVFS1;
PLVBN = PLVFS1;

PLIFS1 = PLVFS1;
PLIFS2 = PLVFS1;
PLIFN = PLVFS1;
PLIBS1 = PLVFS1;
PLIBS2 = PLVFS1;
PLIN = PLVFS1;

WPLIFS1 = PLVFS1;
WPLIFS2 = PLVFS1;
WPLIFN = PLVFS1;
WPLIBS1 = PLVFS1;
WPLIBS2 = PLVFS1;
WPLIN = PLVFS1;

CFS1 = PLVFS1;
CFS2 = PLVFS1;
CFN = PLVFS1;
CBS1 = PLVFS1;
CBS2 = PLVFS1;
CN = PLVFS1;

% Take Hilbert Transform of recordings ERP.erp
% Slice ERP.erp [Channels x time x trial] into a new var [time x channel]
for i = 1:size(ERP.erp,3)
    triali = ERP.erp(:,:,i);
    trialiHT = hilbert(triali);
    [phi, phidiff, phipair, PLIpair, WPLIpair,Cpair] = calcSync(trialiHT,CH);
    % Better name would be PLVtable but phit is phat!
    
    phit(:,:,i) = phipair;
    PLIt(:,:,i) = PLIpair;
    WPLIt(:,:,i) = WPLIpair;
    Ct(:,:,i) = Cpair;
end

% Average over trials
% Manual logical indexing using for loop
rowsS = size(SCORE.F_S1,1);
colsS = size(SCORE.F_S1,2);

%itrial = 1;
% for i = 1:rowsS
%     for j = 1:colsS
%         trial = i*(rowsS-1) + j
%         if (SCORE.F_S1(i,j)==1)
%             phitFS1 = cat(3,phitFS1, phit(:,:,trial));

% % SUPER MEANINGLESS CRAP
 xFS1 = nanmean(nanmean(ERP.erp(:,:,SCORE.F_S1(:)),3),1);
% xFS2 = nanmean(nanmean(ERP.erp(:,:,SCORE.F_S2(:)),3),1);
% xFN = nanmean(nanmean(ERP.erp(:,:,SCORE.F_N(:)),3),1);
% xBS1 = nanmean(nanmean(ERP.erp(:,:,SCORE.B_S1(:)),3),1);
% xBS2 = nanmean(nanmean(ERP.erp(:,:,SCORE.B_S2(:)),3),1);
% xBN = nanmean(nanmean(ERP.erp(:,:,SCORE.B_N(:)),3),1);
% 
 dxFS1 = nanstd(nanstd(ERP.erp(:,:,SCORE.F_S1(:)),0,3),0,1);
% dxFS2 = nanstd(nanstd(ERP.erp(:,:,SCORE.F_S2(:)),0,3),0,1);
% dxFN = nanstd(nanstd(ERP.erp(:,:,SCORE.F_N(:)),0,3),0,1);
% dxBS1 = nanstd(nanstd(ERP.erp(:,:,SCORE.B_S1(:)),0,3),0,1);
% dxBS2 = nanstd(nanstd(ERP.erp(:,:,SCORE.B_S2(:)),0,3),0,1);
% dxBN = nanstd(nanstd(ERP.erp(:,:,SCORE.B_N(:)),0,3),0,1);

 xuFS1 = xFS1+dxFS1;
 xdFS1 = xFS1-dxFS1;
% 
% xuFS2 = xFS2+dxFS2;
% xdFS2 = xFS2-dxFS2;
% 
% xuFN = xFN+dxFN;
% xdFN = xFN-dxFN;
% 
% xuBS1 = xBS1+dxBS1;
% xdBS1 = xBS1-dxBS1;
% 
% xuBS2 = xBS2+dxBS2;
% xdBS2 = xBS2-dxBS2;
% 
% xuBN = xBN+dxBN;
% xdBN = xBN-dxBN;

phitFS1 = phit(:,:,SCORE.F_S1(:));
phitFS2 = phit(:,:,SCORE.F_S2(:));
phitFN = phit(:,:,SCORE.F_N(:));
phitBS1 = phit(:,:,SCORE.B_S1(:));
phitBS2 = phit(:,:,SCORE.B_S2(:));
phitBN = phit(:,:,SCORE.B_N(:));

PLIFS1 = PLIt(:,:,SCORE.F_S1(:));
PLIFS2 = PLIt(:,:,SCORE.F_S2(:));
PLIFN = PLIt(:,:,SCORE.F_N(:));
PLIBS1 = PLIt(:,:,SCORE.B_S1(:));
PLIBS2 = PLIt(:,:,SCORE.B_S2(:));
PLIBN = PLIt(:,:,SCORE.B_N(:));

WPLIFS1 = WPLIt(:,:,SCORE.F_S1(:));
WPLIFS2 = WPLIt(:,:,SCORE.F_S2(:));
WPLIFN = WPLIt(:,:,SCORE.F_N(:));
WPLIBS1 = WPLIt(:,:,SCORE.B_S1(:));
WPLIBS2 = WPLIt(:,:,SCORE.B_S2(:));
WPLIBN = WPLIt(:,:,SCORE.B_N(:));

CFS1 = Ct(:,:,SCORE.F_S1(:));
CFS2 = Ct(:,:,SCORE.F_S2(:));
CFN = Ct(:,:,SCORE.F_N(:));
CBS1 = Ct(:,:,SCORE.B_S1(:));
CBS2 = Ct(:,:,SCORE.B_S2(:));
CBN = Ct(:,:,SCORE.B_N(:));

phitFS1m = nanmean(phitFS1,3);
phitFS2m = nanmean(phitFS2,3);
phitFNm = nanmean(phitFN,3);
phitBS1m = nanmean(phitBS1,3);
phitBS2m = nanmean(phitBS2,3);
phitBNm = nanmean(phitBN,3);

PLIFS1m = nanmean(PLIFS1,3);
PLIFS2m = nanmean(PLIFS2,3);
PLIFNm = nanmean(PLIFN,3);
PLIBS1m = nanmean(PLIBS1,3);
PLIBS2m = nanmean(PLIBS2,3);
PLIBNm = nanmean(PLIBN,3);

WPLIFS1m = nanmean(WPLIFS1,3);
WPLIFS2m = nanmean(WPLIFS2,3);
WPLIFNm = nanmean(WPLIFN,3);
WPLIBS1m = nanmean(WPLIBS1,3);
WPLIBS2m = nanmean(WPLIBS2,3);
WPLIBNm = nanmean(WPLIBN,3);

CFS1m = nanmean(Ct(:,:,SCORE.F_S1),3);
CFS2m = nanmean(Ct(:,:,SCORE.F_S2),3);
CFNm = nanmean(Ct(:,:,SCORE.F_N),3);
CBS1m = nanmean(Ct(:,:,SCORE.B_S1),3);
CBS2m = nanmean(Ct(:,:,SCORE.B_S2),3);
CBNm = nanmean(Ct(:,:,SCORE.B_N),3);

names = {'Fp1'; 'AF7'; 'AF3'; 'F1'; 'F3'; 'F5'; 'F7'; 'FT7'; 'FC5'; 'FC3'; 'FC1'; 'C1'; 'C3'; 'C5'; 'T7'; 'TP7'; 'CP5'; 'CP3'; 'CP1'; 'P1'; 'P3'; 'P5'; 'P7'; 'P9'; 'P07'; 'P03'; '01'; 'Iz'; '0z'; 'P0z'; 'Pz'; 'CPz'; 'Fpz'; 'Fp2'; 'AF8'; 'AF4'; 'Afz'; 'Fz'; 'F2'; 'F4'; 'F6'; 'F8'; 'FT8'; 'FC6'; 'FC4'; 'FC2'; 'FCz'; 'Cz'; 'C2'; 'C4'; 'C6'; 'T8'; 'TP8'; 'CP6'; 'CP4'; 'CP2'; 'P2'; 'P4'; 'P6'; 'P8'; 'P10'; 'P08'; 'P04'; '02'; '1'; '2'; '3';};



% Calculate Instantaneous Phase (if it returns wrapped (-pi to pi), unwrap it (-inf to inf))
% Leave it wrapped, radian

% In all metrics, when computing expected value, assume relative phase is
% uniformly distributed on circular distribution. So E[X(phi)] = mean(phi)

% Plot for proof of concept
figure(1);
trial = 4;
chi = 12;
chj = 24;
t = 1:size(ERP.erp,2);
plot(t,squeeze(ERP.erp(chi,:,trial)));
hold on;
plot(t,squeeze(ERP.erp(chj,:,trial)));
xlabel('Time [ms]');
ylabel('Amplitude [mV]');
legend(['Channel: ',num2str(chi)],['Channel: ',num2str(chj)]);
title(['EEG (PLV: ',num2str(phit(chi,chj,trial)),')']);
hold off;

save('Vars.mat');

figure(1);
t = 1:ts;
plot(t,xFS1);
hold on;
fill([t fliplr(t)],[xuFS1 fliplr(xdFS1)],[.9 .9 .9],'linestyle','none');
xlabel('Time [ms]');
ylabel('Amplitude [mV]');
legend('');
title('Focus Supertarget 1');
hold off;

% figure(2);
% plot(t,xFS2);
% hold on;
% fill([t fliplr(t)],[xuFS2 fliplr(xdFS2)],[.9 .9 .9],'linestyle','none');
% xlabel('Time [ms]');
% ylabel('Amplitude [mV]');
% legend('');
% title('Focus Supertarget 2');
% hold off;

figure(2);
%pcolor(phitFS1m);
imagesc(phitFS1m);
colormap;
title('PLV Focus Supertarget 1');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(3);
imagesc(phitFS2m);
colormap;
title('PLV Focus Supertarget 2');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(4);
imagesc(phitFNm);
colormap;
title('PLV Focus No Target');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

% figure(5);
% imagesc(phitBS1m);
% colormap;
% title('PLV Broad Supertarget 1');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(6);
% imagesc(phitBS2m);
% colormap;
% title('PLV Broad Supertarget 2');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(7);
% imagesc(phitBNm);
% colormap;
% title('PLV Broad No Target');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);

figure(8);
imagesc(PLIFS1m);
colormap;
title('PLI Focus Supertarget 1');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(9);
imagesc(PLIFS2m);
colormap;
title('PLI Focus Supertarget 2');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(10);
imagesc(PLIFNm);
colormap;
title('PLI Focus No Target');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

% figure(11);
% imagesc(PLIBS1m);
% colormap;
% title('PLI Broad Supertarget 1');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(12);
% imagesc(PLIBS2m);
% colormap;
% title('PLI Broad Supertarget 2');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(13);
% imagesc(PLIBNm);
% colormap;
% title('PLI Broad No Target');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);

figure(14);
imagesc(WPLIFS1m);
colormap;
title('WPLI Focus Supertarget 1');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(15);
imagesc(WPLIFS2m);
colormap;
title('WPLI Focus Supertarget 2');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

figure(16);
imagesc(WPLIFNm);
colormap;
title('WPLI Focus No Target');
set(gca,'xtick',[1:67],'xticklabel',names);
set(gca,'ytick',[1:67],'yticklabel',names);

% figure(17);
% imagesc(WPLIBS1m);
% colormap;
% title('WPLI Broad Supertarget 1');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(18);
% imagesc(WPLIBS2m);
% colormap;
% title('WPLI Broad Supertarget 2');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);
% 
% figure(19);
% imagesc(WPLIBNm);
% colormap;
% title('WPLI Broad No Target');
% set(gca,'xtick',[1:67],'xticklabel',names);
% set(gca,'ytick',[1:67],'yticklabel',names);