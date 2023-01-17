function [T,ICC] = plot_ANOVA_Munc13(fileName,shiftFactor,jitterFactor)

% plot_ANOVA_Munc13.m: visualize ANOVA structure in data
%
% [T,ICC] = plot_ANOVA_Munc13(fileName,shiftFactor,jitterFactor)
%
% ex. [T,ICC] = plot_ANOVA_Munc13('dataset4- AP evoked EPSC.xlsx',0.2,0.1);
%
% Inputs:
% - fileName: name of a .xlsx file
% - shiftFactor: how much to shift x-values for cells in a batch
% - jitterFactor: how much to jitter x-values
%
% Outputs:
% - T, a cell array containing the 4 ANOVA tables
% - ICC, an array containing intra-cluster correlations
%
% RTB wrote it, 25 October 2022, gray, drizzly morning
% Originally part of the script, non_hBS_Munc13.m
% RTB updated to include calculation of ICC, 27 October 2022

% To plot a nice looking ANOVA table from the cell array, T:
% digits = [-1 -1 -1 -1 2 2 4];
% maintitle = getString(message('stats:anovan:NWayANOVA'));
% header = getString(message('stats:anovan:AnalysisOfVariance'));
% footer = '';
% % ANOVA Table for Group A
% statdisptable(T(:,:,1),maintitle,header,footer,digits);

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

% In this function, we'll do some statistics to examine structure in the
% data introduced by the hierarchical, nested nature of the experiment. We
% start with a very basic 2-way ANOVA, then proceed to measure the
% 'intra-cluster correlation' values, which are obviously related to the
% ANOVA F-statistic.

% ICC calculations inspired by:
% "A solution to dependency: using multilevel analysis to accommodate nested
% data" Emmeke Aarts, Matthijs Verhage, Jesse V Veenvliet, Conor V
% Dolan & Sophie van der Sluis; Nat. Neurosci. 2014, 17(4):491-496.
%
% PDF: aarts-etal-hierarchical_nested_data_statistics-nn2014.pdf

% Variance. Estimate of the variability in a data set. In nested data, the
% total variance (VarT) is the sum of the variance within research objects
% (VarW, variability among observations taken from the same cell) and the
% variance between research objects (VarB, the variation in cell means).

% Intracluster correlation (ICC). Index of the relative similarity of
% observations taken from the same research object (for example, cell), and
% an indicator of the amount of dependence in the data. The ICC is
% calculated as VarB/[VarB + VarW]. Increasing the differences between
% research objects (VarB) and/or decreasing the differences among measures
% within a research object (VarW) increases the ICC. Experimental
% manipulations (for example, genotype) can increase the variability
% between objects (VarB) and thereby increase the ICC. The part of the ICC
% that can be attributed to the experimental manipulation is called the
% explained ICC. The remainder is called the unexplained ICC.

%% Read the Excel file into a table

if nargin < 1, fileName = 'dataset4- AP evoked EPSC.xlsx'; end
if nargin < 2, shiftFactor = 0.2; end  
if nargin < 3, jitterFactor = 0; end

ds = readtable(fileName);
T = cell(5,7,4);    % Cell array to hold the 4 ANOVA tables
ICC = zeros(4,4);   % Array to hold the 16 ICC values
% Check the column names
%varNames = ds.Properties.VariableNames;

%% Visualize the structure: gscatter

% PK wants to shift cells belonging to each batch:
% The idea is to convert the batches [1,2,3] to [-1,0,1] and then multiply
% these by a shift factor. So, for example, if the shift factor were 0.15,
% cell #1 from batch #1 would have an x-value of 0.85; cell #1 from batch
% #2 would be 1.0 and cell #1 from batch #3 would be 1.15.
% Make a copy to do the shifting for plotting ONLY.
ds2 = ds;
ds2.Cell = ds2.Cell + ((ds2.Batch - 2) .* shiftFactor);

nStrains = 2;   % S_1 = RIMs/ELKS, S_2 = RIMs/ELKS/Munc13
nConds = 2;     % C_1 = KO, C_2 = control
allGroups = {'Group A','Group B','Group C','Group D'};
figure('Name',fileName,'Position', [70 225 1400 670]);

for k = 1: nStrains
    for m = 1:nConds
        
        % select experimental group:
        dsGrp = ds((ds.Strain == k) & (ds.Condition == m),:);
        dsGrp2 = ds2((ds2.Strain == k) & (ds2.Condition == m),:);
        grpNum = (k - 1)*2 + m;
        
        % run 2-way ANOVA:
        [~,tbl,~] = anovan(dsGrp.PSC,{dsGrp.Batch,dsGrp.Cell},...
            'varnames',{'Batch','Cell'},'display','off');
        T(:,:,grpNum) = tbl;
        
        % Calculate the intracluster correlation (ICC) w/r/t/ cell. For
        % each experimental group, we'll actually get 4 values: 3 ICCs for
        % cell within each of the 3 batches, then 1 ICC for batch.
        allBatches = unique(dsGrp.Batch);
        nBatches = length(allBatches);
        ICC_cells = zeros(nBatches,1);
        
        % Note that for this calculation, we want to normalize the variance
        % by 'n', not 'n-1' (default), so we have to specify this when we
        % call 'var' by setting the 2nd argument to '1'. But note also that
        % normalizing by 'n-1' does not have a huge effect on ICC.
        varNl = 1;
        
        % ICC for batches
        % Calculate the mean within-group variance
        varSum = 0;
        grpMeans = zeros(nBatches,1);
        for n = 1:nBatches
            y = dsGrp.PSC(dsGrp.Batch == n);
            varSum = varSum + var(y,varNl);
            grpMeans(n) = mean(y);
        end
        varW = varSum / nBatches;
        
        % Calculate the between-group variance:
        varB = var(grpMeans,varNl);
        ICC_batch = varB / (varB + varW);
        
        % Now calculate the ICC of cells within each batch:
        for n = 1:nBatches
            dsBatch = dsGrp(dsGrp.Batch == n,:);
            allCells = unique(dsBatch.Cell);
            nCells = length(allCells);
            varSum = 0;
            grpMeans = zeros(nCells,1);
            for c = 1:nCells
                y = dsBatch.PSC(dsBatch.Cell == c);
                varSum = varSum + var(y,varNl);
                grpMeans(c) = mean(y);
            end
            varW = varSum / nCells;
            varB = var(grpMeans,varNl);
            ICC_cells(n) = varB / (varB + varW);
        end
        
        % Assign ICC values to the ICC array for return:
        ICC(grpNum,1) = ICC_batch;
        ICC(grpNum,2:end) = ICC_cells;
        
        % Use the title to indicate significance for batch and/or cell
        if tbl{2,7} < 0.05
            bStr = '*';
        else
            bStr = '^o';
        end
        
        if tbl{3,7} < 0.05
            cStr = '*';
        else
            cStr = '^o';
        end
        tStr = [allGroups{grpNum} bStr cStr '  ICC_B:' num2str(ICC_batch,2) ...
            ', ICC_{C1}:' num2str(ICC_cells(1),2) ...
            ', ICC_{C2}:' num2str(ICC_cells(2),2) ...
            ', ICC_{C3}:' num2str(ICC_cells(3),2)];
        
        subplot(2,2,grpNum);
        
        % For some reason, gscatter won't work using variables directly
        % from the table, so we need to break them out:
        x = dsGrp2.Cell;
        y = dsGrp2.PSC;
        g = dsGrp2.Batch;
        color = lines(6); % Generate color values
        
        % jitter x data values for better visibility
        jX = x + ((rand(size(x)) .* jitterFactor) - jitterFactor/2);
        
        % suppress legend:
        %h = gscatter(jX,y,g,color(4:6,:),'dsv',[],'off');
        % with legend:
        h = gscatter(jX,y,g,color(4:6,:),'dsv');
        for z = 1:3
            set(h(z),'LineWidth',1);
        end
        legend('Culture 1','Culture 2','Culture 3');
        title(tStr);
        xlabel('Cell #');
        ylabel('PSC (pA)'); % guessing at units
        set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
    end
end