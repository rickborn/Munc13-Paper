% non_hBS_Munc13.m: non-hierarchical bootstrap of Kaeser lab data
%
% RTB wrote it, 08 October 2022, MKE post bike ride with Tosa Spokesmen

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
    
%% Use 2-way ANOVA to compare within vs. between group variances

[pA,tblA,statsA] = anovan(dsGrpA.PSC,{dsGrpA.Batch,dsGrpA.Cell},...
    'varnames',{'Batch','Cell'});

[pB,tblB,statsB] = anovan(dsGrpB.PSC,{dsGrpB.Batch,dsGrpB.Cell},...
    'varnames',{'Batch','Cell'});

[pC,tblC,statsC] = anovan(dsGrpC.PSC,{dsGrpC.Batch,dsGrpC.Cell},...
    'varnames',{'Batch','Cell'});

[pD,tblD,statsD] = anovan(dsGrpD.PSC,{dsGrpD.Batch,dsGrpD.Cell},...
    'varnames',{'Batch','Cell'});

%% Try to visualize the structure: gscatter

nStrains = 2;   % S_1 = RIMs/ELKS, S_2 = RIMs/ELKS/Munc13
nConds = 2;     % C_1 = KO, C_2 = control
allGroups = {'Group A','Group B','Group C','Group D'};
figure('Name',fileName,'Position', [70 225 990 670]);

for k = 1: nStrains
    for m = 1:nConds
        
        % select experimental group:
        dsGrp = ds((ds.Strain == k) & (ds.Condition == m),:);
        grpNum = (k - 1)*2 + m;
        
        % run 2-way ANOVA:
        [~,tbl,~] = anovan(dsGrp.PSC,{dsGrp.Batch,dsGrp.Cell},...
            'varnames',{'Batch','Cell'},'display','off');
        
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
        x = dsGrp.Cell;
        % jitter x
        jFactor = 0.4;
        jX = x + ((rand(size(x)) .* jFactor) - jFactor/2);
        y = dsGrp.PSC;
        g = dsGrp.Batch;
        color = lines(6); % Generate color values
        % suppress legend:
        %h = gscatter(jX,y,g,color(4:6,:),'dsv',[],'off');
        % with legend:
        h = gscatter(jX,y,g,color(4:6,:),'dsv');
        for z = 1:3
            set(h(z),'LineWidth',1);
        end
        legend('Batch 1','Batch 2','Batch 3');
        title(tStr);
        xlabel('Cell #');
        ylabel('PSC (pA)'); % guessing at units
        set(gca,'LineWidth',1,'FontSize',12,'TickDir','out');
    end
end

%% Non-Hierarchical bootstrap

% For each group, we simply pool all of the data, resample and take the
% mean.

tic;
% For nBoot = 100,000, run time of about 
% For nBoot = 10,000 the run time was around 12 seconds
% For nBoot = 1,000, run time is 1.5 seconds

% define constants:
nBoot = 10000;

% variable to store the grand mean value of the EPSCs for each group (rows)
% and each bootstrap iteration (columns):
allMeans = zeros(nStrains*nConds, nBoot);

for thisBoot = 1:nBoot
    for thisStrain = 1: nStrains
        for thisCond = 1:nConds
            
            % Grab the subset of the data corresponding to this group:
            dsGrp = ds((ds.Strain == thisStrain) & (ds.Condition == thisCond),:);
            
            % resample with replacement from all recordings
            nGrp = height(dsGrp);
            dataStar = dsGrp.PSC(unidrnd(nGrp,nGrp,1));
            
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

figure
histogram(Tboot);
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

% make it pretty
xlabel('T^{*}');
ylabel('Frequency');
title([fileName '     p(H0|Data) = ' num2str(pValue,3)]);
legend([h1,h2],{'Experimental value','95% CI'});
set(gca,'FontSize',12);