% perm_test_Munc13.m: permutation test of Kaeser lab data
%
% RTB wrote it, 02 October 2022, Sunday morning prior to Lake Michigan swim

% When the Kaeser lab knocks out two genes involved in synaptic
% release (RIMs & ELKS), synaptic transmission is reduced by about 80%. Now
% they want to ask whether KO'ing a 3rd gene, Munc13, further reduces
% synaptic transmission. However, since the double KO (RIMs/ELKS) is in a
% different mouse strain from the triple KO (RIMs/ELKS/Munc13), the
% synaptic transmission in each strain's KO must be compared to its own
% control.
%
% S_1 = Strain 1: RIMs/ELKS-Cre
% S_2 = Strain 2: RIMs/ELKS/Munc13-Cre
%
% To perform an experiment in a given strain, several mice are killed, the
% hippocampi (HC) are dissected, cells dissociated and pooled together in
% one big primary culture. This culture is then divided into two groups:
% one is treated with a lentivirus containing Cre (= KO); the other is
% treated with a lentivirus containing a dead Cre, 'delta-Cre' (= Control).
% From each culture, multiple cells are tested. Each cell is patch clamped
% and synaptic transmission is tested by measuring the size of the
% post-synaptic current (EPSC) evoked by an action potential. This
% measurment is generally repeated 5 times for each cell.
%
% C_1 = Condition 1: KO
% C_2 = Condition 2: Control
%
% Thus each group of 5 related measurements is uniquely identified by FOUR
% numbers, which will be variables (columns) in the Excel file. The Excel
% file will be read in as a Table variable, 'ds':
%
% ds.Strain [1:2]
% ds.Cond [1:2]
% ds.Batch [1:nBatches], may vary for different experiments
% ds.Cell [1:nCells], may vary for different batches
%
% For ease of discussion, we'll think of four groups:
%
% Group A: S_1, C_1 (RIMS/ELKS KO)
% Group B: S_1, C_2 (RIMS/ELKS Control)
% Group C: S_2, C_1 (RIMS/ELKS/Munc13 KO)
% Group D: S_2, C_2 (RIMS/ELKS/Munc13 Control)
%
% The scientific question is whether knocking out Munc13 (triple mutant,
% S2) causes a greater relative decrease in synaptic transmission as
% compared to the double mutant (S1). So we will define our test statistic,
% T, as:
%
% T = (mean(Group_A) / mean(Group_B)) / . . . 
%     (mean(Group_C) / mean(Group_D)) 
%
% The null value for this statistic is 1; our alternate hypothesis is:
%       T > 1

% This code was adapted from 'hBS_Munc13.m' in which I performed a
% hierarchical bootstrap on the test statistic, T. Here, we will do a
% simpler permutation test to simulate the behavior of T under H0. To do
% this, I'll simply scramble the labels corresponding to 'Strain' and
% repeat the calculation of T_star.

%% Read the Excel file into a table

%fileName ='fakedata1.xlsx';
fileName = 'dataset 1-sucrose IPSC-v4.xlsx';
%fileName = 'dataset 2-sucrose E-v1.xlsx';
%fileName = 'dataset 3- AP evoked IPSC.xlsx';
ds = readtable(fileName);

% Check the column names
varNames = ds.Properties.VariableNames;

%% Calculate the actual value of our test statistic, T

dsGrpA = ds((ds.Strain == 1) & (ds.Condition == 1),:);
dsGrpB = ds((ds.Strain == 1) & (ds.Condition == 2),:);
dsGrpC = ds((ds.Strain == 2) & (ds.Condition == 1),:);
dsGrpD = ds((ds.Strain == 2) & (ds.Condition == 2),:);

Treal = (mean(dsGrpA.PSC,'omitnan') / mean(dsGrpB.PSC,'omitnan')) / ...
        (mean(dsGrpC.PSC,'omitnan') / mean(dsGrpD.PSC,'omitnan'));

%% Permutation Test

% # of observations for each strain
nStrain1 = sum(ds.Strain == 1);
nStrain2 = sum(ds.Strain == 2);
nTotal = nStrain1 + nStrain2;

% time it
tic;
% This test is much faster:
% For 100,000 permutations, it took just 2.5 min.
% define constants:
nPerm = 10000;

% variable to store the grand mean value of the EPSCs for each group (rows)
% and each bootstrap iteration (columns):
Tperm = zeros(nPerm,1);

for k = 1:nPerm
    % Shuffle the strain labels:
    ds.Strain = ds.Strain(randperm(nTotal));
        
    % redo the calculation:
    dsGrpA = ds((ds.Strain == 1) & (ds.Condition == 1),:);
    dsGrpB = ds((ds.Strain == 1) & (ds.Condition == 2),:);
    dsGrpC = ds((ds.Strain == 2) & (ds.Condition == 1),:);
    dsGrpD = ds((ds.Strain == 2) & (ds.Condition == 2),:);
    
    Tperm(k,1) = (mean(dsGrpA.PSC,'omitnan') / mean(dsGrpB.PSC,'omitnan')) / ...
                 (mean(dsGrpC.PSC,'omitnan') / mean(dsGrpD.PSC,'omitnan'));
end

elapsedTimeInSeconds = toc;

%% Calculate p-values

% calculate 1-tailed p-value:
pValue1t = sum(Tperm >= Treal) / nPerm;

if pValue1t == 0
    pValue1t = 1 / (nPerm + 1);
end

% calculate 2-tailed p-value:
pValue2t = sum((Tperm >= Treal) | (Tperm <= (1/Treal))) / nPerm;

if pValue2t == 0
    pValue2t = 1 / (nPerm + 1);
end

%% Plot a histogram of our H0 distribution

%figure
histogram(Tperm);
hold on

% add a line for our experimental value
ax = axis;
h1 = line([Treal, Treal],[ax(3), ax(4)]);
set(h1, 'Color', [0.6350, 0.0780, 0.1840],'LineWidth',2);
h2 = line([1/Treal, 1/Treal],[ax(3), ax(4)]);
set(h2, 'Color', [0.6350, 0.0780, 0.1840],'LineWidth',2);

% make it pretty
xlabel('T_{perm}');
ylabel('Frequency');
title([fileName '     p-value = ' num2str(pValue2t,3)]);
legend(h1,{'Experimental value'});
set(gca,'FontSize',12);