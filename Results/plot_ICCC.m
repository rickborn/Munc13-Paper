% plot_ICC.m: plot histogram of intra-cluster correlations
%
% RTB wrote it, 31 October 2022, listening to Supreme Court arguments on
% SFFA vs. Harvard (community Zoom meeting)

% load data
% The values in all_ICC.mat were returned by 'plot_ANOVA_Munc13.m'
load all_ICC.mat

% Data in array named 'ICC'
% For each of the 6 data files (3rd dimension), the 1st column contains the
% ICCs for 'batch' and the 2nd through 4th columns contain the ICCs for the
% 'cells' in each 'batch'. Experiments 1 and 2 (sucrose) only have one
% measurement per cell, so they have (meaningless) values of ICC_cell of 1.
% We should omit these from the cell values.

% separate out the ICCs for 'batch' and for 'cell'
all_ICC_Batch = ICC(:,1,:);
all_ICC_Cells = ICC(:,2:end,3:end);

% plot histograms:
hB = histogram(all_ICC_Batch(:));
hold on
hC = histogram(all_ICC_Cells(:));

% You'll need to adjust the bin sizes using:
% nBins = morebins(hC)
% nBins = fewerbins(hC)

% Label axes:
xlabel('Intra-cluster correlation');
ylabel('Frequency');

set(gca,'FontSize',12,'TickDir','out');
legend({'Batch','Cell'},'Location','north');