% hBS_Munc13_noBatch.m: hBS of Kaeser lab data, omitting batch
%
% RTB wrote it, 29 October 2022, after walk at Millenium with Anne & Nena
% Adapted from hBS_Munc13.m
%
% In this version, we'll pool all of the cells from a given batch and
% resample from cells and sweeps only. The motivation for this is:
%   1. Resampling from n = 3 batches is probably not legit
%   2. The ICC for batches tends to be low.

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

%% Read the Excel file into a table

%fileName ='fakedata1.xlsx';
%fileName = 'dataset 1-sucrose IPSC-v4.xlsx';
%fileName = 'dataset 2-sucrose E-v1.xlsx';
%fileName = 'dataset 3- AP evoked IPSC.xlsx';
fileName = 'dataset4- AP evoked EPSC.xlsx';
%fileName = 'dataset4- AP evoked EPSC withNaNs.xlsx';
%fileName = 'dataset 5-EPSC-PPR(50ms).xlsx';
%fileName = 'dataset 6-IPSC-PPR(20ms).xlsx';
ds = readtable(fileName);

% Check the column names
varNames = ds.Properties.VariableNames;

%% Calculate the actual value of our test statistic, T

dsGrpA = ds((ds.Strain == 1) & (ds.Condition == 1),:);  % double KO Cre
dsGrpB = ds((ds.Strain == 1) & (ds.Condition == 2),:);  % double KO control
dsGrpC = ds((ds.Strain == 2) & (ds.Condition == 1),:);  % triple KO Cre
dsGrpD = ds((ds.Strain == 2) & (ds.Condition == 2),:);  % triple KO control

Treal = (mean(dsGrpA.PSC,'omitnan') / mean(dsGrpB.PSC,'omitnan')) / ...
        (mean(dsGrpC.PSC,'omitnan') / mean(dsGrpD.PSC,'omitnan'));

%% Re-code cells

% Idea is to give each cell within a given experimental condition a unique
% cell #. Each cell's ID is thus defined by two numbers: batch and cell.
allCells = unique([dsGrp.Batch,dsGrp.Cell],'rows');

% Create a new column for the unique cell #:
dsGrp.UCell = zeros(height(dsGrp),1);

% Loop through allCells and assign unique cell #.
cellID = 1;
for k = 1:length(allCells)
    dsGrp.UCell((dsGrp.Batch == allCells(k,1)) & (dsGrp.Cell == allCells(k,2))) = cellID;
    cellID = cellID + 1;
end

% Note: Doing it like this will create gaps in unique cell #, because
% different batches can have different #s of cells. Can solve this by doing
% the recoding withing experimental group.

%% Hierarchical bootstrap

% For each group, we first resample (with replacement) from the batches.
% Then, for each of these batches, we resample from the cells; then for
% each cell we resample from the (usually 5) technical replicates. For the
% "sucrose" experiments, there is only one data point per cell.

% Many 'for' loops, and a lot of subsetting of the data (which can probably
% be avoided with more thought), so this might be kinda slow. Let's time it
% for different values of nBoot
tic;
% For nBoot = 100,000, run time of about 45 min.
% For nBoot = 10,000 the run time was around 353 seconds (almost 6 minutes)
% For nBoot = 1,000, run time is 28 seconds
% For nBoot = 100, run time is about 3 seconds

% define constants:
nBoot = 10;
nStrains = 2;   % S_1 = RIMs/ELKS, S_2 = RIMs/ELKS/Munc13
nConds = 2;     % C_1 = KO, C_2 = control

% variable to store the grand mean value of the EPSCs for each group (rows)
% and each bootstrap iteration (columns):
allMeans = zeros(nStrains*nConds, nBoot);

for thisBoot = 1:nBoot
    for thisStrain = 1: nStrains
        for thisCond = 1:nConds
            
            % Temporary variable to hold the resampled EPSC values:
            dataStar = [];
            
            % Grab the subset of the data corresponding to this group:
            dsGrp = ds((ds.Strain == thisStrain) & (ds.Condition == thisCond),:);
            
            % How many batches for this group?
            nBatches = length(unique(dsGrp.Batch,'rows'));
            
            % re-sample with replacement from batches (Note that this
            % assumes batches go from 1:nBatches with no gaps.)
            bStarNdx = unidrnd(nBatches, nBatches, 1);
            
            for thisBatch = 1:nBatches
                dsBatch = dsGrp(dsGrp.Batch == bStarNdx(thisBatch),:);
                nCells = length(unique(dsBatch.Cell,'rows'));
                cStarNdx = unidrnd(nCells, nCells, 1);
                
                for thisCell = 1:nCells
                    % dsCell contains all the measurements for one cell:
                    dsCell = dsBatch(dsBatch.Cell == cStarNdx(thisCell),:);
                    %nSweeps = length(unique(dsCell.PSC,'rows'));
                    nSweeps = height(dsCell);
                    
                    if nSweeps == 1
                        dataStar = [dataStar; dsCell.PSC];
                    else
                        swStarNdx = unidrnd(nSweeps, nSweeps, 1);
                        dataStar = [dataStar; dsCell.PSC(swStarNdx)];
                    end
                end
            end
            
            % Store mean of dataStar in allMeans
            allMeans(((thisStrain - 1)*2) + thisCond, thisBoot) = mean(dataStar,'omitnan');
            
        end
    end
end


% Calculate our test statistic, T, for each bootstrap iteration:
Tboot = (allMeans(1,:) ./ allMeans(2,:)) ./ (allMeans(3,:) ./ allMeans(4,:));

elapsedTimeInSeconds = toc;

%% Calculate confidence intervals and a p-value

% calculate 95% CI using percentile method:
sortedTboot = sort(Tboot);
[idxHi,idxLo] = cindx(0.05,nBoot);
CIhi = sortedTboot(idxHi);
CIlo = sortedTboot(idxLo);

% calculate p-value:
pValue = sum(Tboot <= 1) / nBoot;

if pValue == 0
    pValue = 1 / (nBoot + 1);
end

%% Plot a histogram of our bootstrap distribution

%figure
h = histogram(Tboot);
hold on

% add a line for our experimental value
ax = axis;
h1 = line([Treal, Treal],[ax(3), ax(4)]);
set(h1, 'Color', [0.6350, 0.0780, 0.1840],'LineWidth',2);

% add dashed lines for confidence intervals
h2 = line([CIhi, CIhi],[ax(3), ax(4)]);
% green lines:
set(h2, 'Color', [0.4660, 0.6740, 0.1880],'LineStyle','--','LineWidth',1.5);
h3 = line([CIlo, CIlo],[ax(3), ax(4)]);
set(h3, 'Color', [0.4660, 0.6740, 0.1880],'LineStyle','--','LineWidth',1.5);

% make it pretty
xlabel('T^{*}');
ylabel('Frequency');
%title([fileName '     p(H0|Data) = ' num2str(pValue,3)]);
legend([h1,h2],{'Experimental value','95% CI'});
set(gca,'FontSize',12);

% For comparison purposes, it's useful to adjust the # of bins:
% nBins = morebins(h)
% nBins = fewerbins(h)