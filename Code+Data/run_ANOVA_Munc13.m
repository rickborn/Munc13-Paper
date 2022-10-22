% run_ANOVA_Munc13.m: wrapper function to analyze multiple data sets
%
% RTB wrote it, 10 October 2022, after spectacular moon-set/sun-rise swim
% at Klode Beach in Whitefish Bay (MKE)

%% Set up storage

% Load in cell array with all of the file names:
load all_data_files_names.mat

% Create a table to store results of analyses:
varNames = {'file_name','grp_name','batch_df','batch_F','batch_pVal',...
    'cell_df','cell_F','cell_pVal'};
varTypes = {'string','string','uint8','double','double',...
    'uint8','double','double'};
grpNames = {'GrpA','GrpB','GrpC','GrpD'};
nStrains = 2;
nConds = 2;
nGroups = length(grpNames);
nFiles = length(allFileNames);
nRows = nFiles * nGroups;
nCols = length(varNames);
color = lines(6); % Generate color values for color grouping of batches

T = table('Size',[nRows,nCols],'VariableTypes',varTypes,'VariableNames',varNames);

%% Run analyses

for k = 1:nFiles
    ds = readtable(allFileNames{k});
    
    if k > 2
        figure('Name',allFileNames{k},'Position', [70 225 990 670]);
    end
    
    for m = 1:nStrains
        for n = 1:nConds
            
            % Select the group
            dsGrp = ds((ds.Strain == m) & (ds.Condition == n),:);
            [~,tbl] = anovan(dsGrp.PSC,{dsGrp.Batch,dsGrp.Cell},...
                'varnames',{'Batch','Cell'},'display','off');
            
            % fancy math to figure out row indexing:
            thisGrp = ((m - 1)*2) + n;
            thisRow = ((k - 1)*4) + thisGrp;
            
            % Assign values to table
            T.file_name(thisRow) = allFileNames{k};
            T.grp_name(thisRow) = grpNames{thisGrp};
            T.batch_df(thisRow) = tbl{2,3};
            T.batch_F(thisRow) = tbl{2,6};
            T.batch_pVal(thisRow) = tbl{2,7};
            T.cell_df(thisRow) = tbl{3,3};
            T.cell_F(thisRow) = tbl{3,6};
            T.cell_pVal(thisRow) = tbl{3,7};
            
            % Plot each experimental group (but not for sucrose expts):
            % Use the title to indicate significance for batch and/or cell
            if k > 2
                if T.batch_pVal(thisRow) < 0.05
                    bStr = '*';
                else
                    bStr = '^o';
                end
                
                if T.cell_pVal(thisRow) < 0.05
                    cStr = '*';
                else
                    cStr = '^o';
                end
                tStr = [grpNames{thisGrp} bStr cStr];
                
                subplot(2,2,thisGrp);
                % For some reason, gscatter won't work using variables directly
                % from the table, so we need to break them out:
                x = dsGrp.Cell;
                % jitter x
                jFactor = 0.3;
                jX = x + ((rand(size(x)) .* jFactor) - jFactor/2);
                y = dsGrp.PSC;
                g = dsGrp.Batch;
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
    end
end


%% Save the results table

resultFileName = 'results_ANOVA';
save(resultFileName, 'T');
    
% Also write to an excel spreadsheet
writetable(T,'results_ANOVA.xlsx');