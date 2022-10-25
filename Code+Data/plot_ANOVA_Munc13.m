function [T] = plot_ANOVA_Munc13(fileName,shiftFactor,jitterFactor)

% plot_ANOVA_Munc13.m: visualize ANOVA structure in data
%
% RTB wrote it, 25 October 2022, gray, drizzly morning
% Originally part of the script, non_hBS_Munc13.m
%
% [T] = plot_ANOVA_Munc13(fileName,shiftFactor,jitterFactor)
%
% Inputs:
% - fileName: name of a .xlsx file
% - shiftFactor: how much to shift x-values for cells in a batch
% - jitterFactor: how much to jitter x-values
%
% Outputs:
% - T, a cell array containing the 4 ANOVA tables

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

% In this script, we'll explore whether doing the bootstrap hierarchically
% matters by running it non-hierarchically and comparing with the results
% of the hierarchical bootstrap. We'll also do some additional statistics
% to examine structure in the data introduced by the nested nature of the
% experiment.

%% Read the Excel file into a table

if nargin < 1, fileName = 'dataset4- AP evoked EPSC.xlsx'; end
if nargin < 2, shiftFactor = 0.2; end  
if nargin < 3, jitterFactor = 0; end

ds = readtable(fileName);
T = cell(5,7,4);    % Cell array to hold the 4 ANOVA tables
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
        tStr = [allGroups{grpNum} bStr cStr];
        
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