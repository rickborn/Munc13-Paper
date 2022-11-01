% plot_Many_Histograms.m: superimpose three related histograms
%
% RTB wrote it, 31 October 2022, listening to KC "The ConstruKction of
% Light"

cd 'C:\usr\rick\doc\Students\Kaeser, Pascal\Munc13 Paper\Munc13-Paper\Results'

load 'dataset4- AP evoked EPSC_hFlag_1_nBoot_100000.mat'
Tb1 = Tb;
load 'dataset4- AP evoked EPSC_noBatch_nBoot_100000.mat'
Tb2 = Tb;
load 'dataset4- AP evoked EPSC_hFlag_0_nBoot_100000.mat'
Tb3 = Tb;

h1 = plot_Tboot_Histogram(Tb1, Tr, 0.05, 6);
h2 = plot_Tboot_Histogram(Tb2, Tr, 0.05, 3);
h3 = plot_Tboot_Histogram(Tb3, Tr, 0.05, 4);

legend 'off'
title('dataset4- AP evoked EPSC.xlsx');
axis([0.5, 3.5, 0, 3500]);