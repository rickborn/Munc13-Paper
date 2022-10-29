% run_all_noBatch_Munc13.m: wrapper function to analyze multiple data sets
%
% This is for the "noBatch" version, so we won't do the non-hBS
%
% RTB wrote it, 29 October 2022

%% Specify analysis parameters

nBoot = 100000;
myAlpha = 0.05;

%% Set up storage

% Load in cell array with all of the file names:
load all_data_files_names.mat

% Create a table to store results of analyses:
varNames = {'file_name','T_real','SE','pVal','CI_lo','CI_hi'};
varTypes = {'string','double','double','double','double','double'};
% We run the bootstrap only once for each file:
nRows = length(allFileNames);
nCols = length(varNames);

T = table('Size',[nRows,nCols],'VariableTypes',varTypes,'VariableNames',varNames);

%% Run analyses

for k = 1:nRows
    
    % Get the relevant file name:
    thisFileName = allFileNames{k};
    
    % Run the analysis
    [Tr,Tb,allMu,SE,CI,p] = hBS_Munc13_noBatch_function(thisFileName,nBoot,myAlpha,0);
    
    % NOTE: If we save T_boot for each analysis, we can reproduce the
    % histograms at our leisure. Consider doing this.
    [pathstr, fileName, fileExt] = fileparts(thisFileName);
    BS_fName = [fileName,'_noBatch_nBoot_',num2str(nBoot),'.mat'];
    save(BS_fName,'Tr','Tb','allMu');
    
    % Assign values to table
    T.file_name(k) = thisFileName;
    T.T_real(k) = Tr;       % experimental value of test statistic, T
    T.SE(k) = SE;           % standard error of T
    T.pVal(k) = p;          % p(H0|data)
    T.CI_lo(k) = CI(1);     % lower confidence interval
    T.CI_hi(k) = CI(2);     % upper confidence interval
end

%% Save the results table

resultFileName = ['results_noBatch_nBoot_' num2str(nBoot) '.mat'];
save(resultFileName, 'T');

% Also write to an excel spreadsheet
excelFileName = ['results_noBatch_nBoot_' num2str(nBoot) '.xlsx'];
writetable(T,excelFileName);