% write_Munc13_results.m: appends mat files and writes to Excel
%
% RTB wrote it, 15 October 2022, chilly Saturday morning

%% Make a list of the results files

d = 'C:\usr\rick\doc\Students\Kaeser, Pascal\Munc13 Paper\Results';
files = dir(fullfile(d, '*.mat'));
% Double check to make sure these are the correct files!
files = files(1:12,1);

% Get the file names out of the struct and into a cell array:
F = struct2cell(files);
fileNames = F(1,:);

%% Create a table to store concatenated data

% Column names will be the name of the results file
varTypes = {'double','double','double','double','double','double',...
    'double','double','double','double','double','double'};
nRows = 100000;
nCols = length(fileNames);

T = table('Size',[nRows,nCols],'VariableTypes',varTypes,'VariableNames',fileNames);

%% Read through list of files, load and paste into correct column

for k = 1:nCols
    load(fileNames{k});
    T_boot = table(Tb');
    T(:,k) = T_boot;
end

%% Table reality check: compare with original p-values:

p = sum(T{:,1} <= 1) / 100000;

%% Write table to an Excel file

writetable(T,'all_hBootstrap_replicates.xlsx');