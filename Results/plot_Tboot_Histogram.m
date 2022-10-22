function [h] = plot_Tboot_Histogram(Tboot, Treal, myAlpha)

% plot_Tboot_Histogram.m: plot histogram of sampling distribution, T*
%
% [h] = plot_Tboot_Histogram(Tboot, myAlpha)
%
% ex. h = plot_Tboot_Histogram(Tboot, 0.05)
%
% Inputs:
% - Tboot: bootstrap replicates, T*
% - Treal: experimental value of T
% - myAlpha, significance level for confidence intervals (default = 0.05)
%
% Outputs:
% - h: handle to histogram
%
% RTB wrote it, 14 October 2022, watching game #2 of NY Yankees vs.
% Cleveland Guardians (ALDS)

%% Calculate standard error, confidence intervals and a p-value

% find # of bootstrap replicates:
nBoot = length(Tboot);

% calculate the standard error:
SE = std(Tboot);

% calculate 95% CI using percentile method:
sortedTboot = sort(Tboot);
idxHi = ceil(nBoot * (1 - myAlpha/2));
idxLo = floor(nBoot * (myAlpha/2));
CIhi = sortedTboot(idxHi);
CIlo = sortedTboot(idxLo);
CI = [CIlo, CIhi];

% calculate p-value:
pValue = sum(Tboot <= 1) / nBoot;

if pValue == 0
    pValue = 1 / (nBoot + 1);
end

%% Plot a histogram of our bootstrap distribution


h = histogram(Tboot);
% h = hist(Tboot,50);   % seems to be better for sharing vectorized
% graphics
hold on

% add a line for our experimental value
ax = axis;
h1 = line([Treal, Treal],[ax(3), ax(4)]);
set(h1, 'Color', [0.6350, 0.0780, 0.1840],'LineWidth',2);

% add dashed lines for confidence intervals
h2 = line([CIhi, CIhi],[ax(3), ax(4)]);
set(h2, 'Color', [0.6350, 0.0780, 0.1840],'LineStyle','--','LineWidth',1.5);
h3 = line([CIlo, CIlo],[ax(3), ax(4)]);
set(h3, 'Color', [0.6350, 0.0780, 0.1840],'LineStyle','--','LineWidth',1.5);

% label axes, etc.
xlabel('T^{*}');
ylabel('Frequency');
title(['p(H0|Data) = ' num2str(pValue,3)]);
ciStr = [num2str((1-myAlpha)*100) '% CI'];
legend([h1,h2],{'Experimental value',ciStr});
set(gca,'FontSize',12,'TickDir','out');
