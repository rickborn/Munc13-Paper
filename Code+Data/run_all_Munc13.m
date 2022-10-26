% run_all_Munc13.m: wrapper function to analyze multiple data sets
%
% RTB wrote it, 09 October 2022, while watching Packers v. Giants (played
% in London)

%% Specify analysis parameters

nBoot = 100000;
myAlpha = 0.05;

%% Set up storage

% Load in cell array with all of the file names:
load all_data_files_names.mat

% Create a table to store results of analyses:
varNames = {'file_name','h_Flag','T_real','SE','pVal','CI_lo','CI_hi'};
varTypes = {'string','logical','double','double','double','double','double'};
% We run the bootstrap twice for each file:
% once hierarchical (odd # rows), once non-hierarchical (even # rows).
nRows = length(allFileNames) * 2;
nCols = length(varNames);

T = table('Size',[nRows,nCols],'VariableTypes',varTypes,'VariableNames',varNames);

%% Run analyses

for k = 1:nRows
    
    % Results of hierarchical bootstrap on odd # rows; non-hBS on even
    hFlag = mod(k,2);
    
    % Get the relevant file name:
    thisFileName = allFileNames{ceil(k/2)};
    
    % Run the analysis
    [Tr,Tb,allMu,SE,CI,p] = hBS_Munc13_function(thisFileName,hFlag,nBoot,myAlpha,0);
    
    % NOTE: If we save T_boot for each analysis, we can reproduce the
    % histograms at our leisure. Consider doing this.
    [pathstr, fileName, fileExt] = fileparts(thisFileName);
    BS_fName = [fileName,'_hFlag_',num2str(hFlag),'_nBoot_',num2str(nBoot),'.mat'];
    save(BS_fName,'Tb','allMu');
    
    % Assign values to table
    T.file_name(k) = thisFileName;
    T.h_Flag(k) = hFlag;    % 1 = hierarchical bootstrap
    T.T_real(k) = Tr;       % experimental value of test statistic, T
    T.SE(k) = SE;           % standard error of T
    T.pVal(k) = p;          % p(H0|data)
    T.CI_lo(k) = CI(1);     % lower confidence interval
    T.CI_hi(k) = CI(2);     % upper confidence interval
end

%% Save the results table

resultFileName = ['results_nBoot_' num2str(nBoot) '.mat'];
save(resultFileName, 'T');

% Also write to an excel spreadsheet
excelFileName = ['results_nBoot_' num2str(nBoot) '.xlsx'];
writetable(T,excelFileName);