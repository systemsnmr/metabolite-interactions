function main(curdir,fsamode,ploton)

%%% Help / annotation
%====================================================
% Example call: main('.\',1,0)
% This function collects info about nmr spectra and detects the peaks from
% peaklist
% 
% Input:
% curdir: working directory, if it is current folder say  '.\'
% fsamode:
%   0 : reference to M spectrum(PMP-M)/M)
%   1 : use old way to determine fsa 
%   2 : compute relaxation factor instead of FSA [default]
% ploton: 1: plot every single peak, 0: don't
%
% Output:
% hits are saved in folder 'results_hits' (or in other folder, as specified)
% figures are saved in 'results_figures' (or in other folder, as specified)

%%% Params for testing (are set only when function is run w/o any parameters
%====================================================
if nargin == 0
    fprintf('No input args provided. Running test example.\n');
    param.curdir = '.\';
    param.fsamode = 2;
    param.ploton = 0;
else
    % translate input arguments to structure
    param.curdir = curdir;
    param.fsamode = fsamode;
    param.ploton = ploton;
    clearvars curdir datasets_file frompeaklist recopy_spectra fsamode ploton;
end

% define additional parameters
% output data
param.out.figures = fullfile(param.curdir,'results\');
param.out.hits = fullfile(param.curdir,'results_matfiles\');
param.out.spectra = fullfile(param.curdir,'input_spectra\');

% input data
param.in.spectra = fullfile(param.curdir,'input_spectra\');
param.in.peaklist = fullfile(param.curdir,'input_peaklists\');
param.in.datasets_file = fullfile(param.curdir,'input_other\datasets.csv');

% load structre with info on spectra names
load(fullfile(param.curdir,'input_other\spectranames.mat'));
param.in.spectranames = spectranames;

% analysis parameters
param.recopy_spectra = 0; % recopy spectra from original destination
param.LocalGlobalThreshold = 5; % decides when to use global or local noise

%%% Error and type-check input arguments
%====================================================
if ~exist(param.in.datasets_file, 'file') == 2
    error('datasets_file has to be a valid filename!');
end
if ~ismember(param.recopy_spectra,[0 1])
    error('recopy_spectra has to be either 0 or 1!');
end
if ~ismember(param.fsamode,[0 1 2])
    error('fsamode has to be either 0, 1 or 2!');
end
if ~ismember(param.ploton,[0 1])
    error('ploton has to be either 0 or 1!');
end


%%% Actual script
%====================================================

% read info about datasets
try
    [~,xout,~] = xlsread(param.in.datasets_file,'A1:E500');
    xout = regexp(xout, ',', 'split');
    xout = vertcat(xout{:});
    raw{1} = xout(:,1);
    raw{2} = xout(:,2);    
catch
    warning('xlsxread didnt work, try textscan instead.');
    fileID = fopen(param.in.datasets_file);
    raw = textscan(fileID,'%s %s','delimiter',',');
    fclose(fileID);
end


% make structure for dataset info
dsets = struct('name', raw{1},'path', cellfun(@(dir,name) fullfile(dir,name), raw{2}, raw{1}, ...
    'UniformOutput', false),'p_names', [], 'm_names', []);

% copy spectra txt if selected
fprintf(1,'Copying processed spectra...\n');
if param.recopy_spectra
    copy_spectra(dsets,param.in.spectra); % make a local copy of spectra for analysis
end

% get more info on datasets
fprintf(1,'Getting protein/sample infos ...\n');
dsets = get_dataset_info(dsets,param.in.spectra); % Get protein/mix names and expnos

% run analysis of all datasets
fprintf(1,'Running analysis...\n');
dsets = run_analysis(dsets,param);

% Tests
for i=1:numel(dsets)
    fprintf(1,'== For dset %s \n', dsets(i).name)
    dsets(i)
end

end

%============================================================
function dsets = run_analysis(dsets,param)
%============================================================
if ~exist(param.out.figures, 'dir')
    mkdir(param.out.figures); % folder for output figures
    disp('Created folder for figure output');
end
if ~exist(param.out.hits, 'dir')
    mkdir(param.out.hits); % folder for hits structures
    disp('Created folder for hits output');
end

n_sets = numel(dsets);

for s = 1:n_sets
    d = dsets(s);

    for p = 1:d.n_prot
        for m = 1:d.n_mix
            
            mix = d.m_names{m};
            protein = d.p_names{p};
            spectraPath = d.spectra;
            
            % peak picking based on peaklist
            main_PeakPicking(spectraPath,mix,protein,param.LocalGlobalThreshold,...
                param.fsamode,param.ploton,param.in, param.out);
            
            close all hidden;
        end        
    end
    
    clear d;
end

end


%============================================================
function copy_spectra(dsets,copypath)
%============================================================

if ~exist(copypath, 'dir')
    mkdir(copypath);
end

n_sets = numel(dsets);

for i=1:n_sets
    
    fprintf(1,'== %s \n', dsets(i).name)
    
    % create directories for each expt
    currdir = fullfile(copypath,dsets(i).name);
    if ~exist(currdir, 'dir')
        mkdir(currdir);
    end
      
    % spectra folder
    spectrapath = fullfile(dsets(i).path,'spectra');  
  
    % extract relevant filenames
    content = dir(spectrapath);
    allfnames = {content.name};
    allfnames = allfnames(strncmp('T1',allfnames,2));
    
    
    % check how many files we have -> if theres one per mix
    if mod(length(allfnames),4) > 0
        disp('Number of filenames is off! expect 1 file per ccm mix');
    end
    
    % dcopy each file
    for j = 1:length(allfnames)
        ff = fullfile(spectrapath,allfnames{j});
        copyfile(ff, currdir);    
    end
        
end

end


%============================================================
function dsets = get_dataset_info(dsets,spectrapath)
%============================================================
% input: dsets structure with path
% output: dsets struct populated with names, expnos, ..

n_sets = numel(dsets);

    
% read content of folder and extract pnames and mixnames from there
for i = 1:n_sets
    
    % get path to txts
    currpath = fullfile(spectrapath,dsets(i).name);
    content = dir(currpath);
    
    % split filenames
    allfnames = {content.name};
    allfnames = allfnames(strncmp('T1',allfnames,2));
    
    % save all filenames
    dsets(i).spectra = currpath;
    dsets(i).allfilenames = allfnames;
   
    
    % get proteins & mixes
    % find all std filenames
    stdnames = {};
    for j = 1:length(allfnames)
        if find(regexpi(allfnames{j},'T1r_\w{3,4}_\w{4}.txt'))
            stdnames(end+1) = allfnames(j);
        end
    end
    splitcont = regexp(stdnames, '_', 'split');
    splitcont = vertcat(splitcont{:});
    proteins = unique(splitcont(:,2));
    mixes = unique(splitcont(:,3));
    
    % remove .txt
    for j = 1:length(mixes)
        mixes{j} = mixes{j}(1:end-4);
    end
      
    dsets(i).n_prot = numel(proteins);
    dsets(i).p_names = proteins';
    
    dsets(i).n_mix = numel(mixes);
    dsets(i).m_names = mixes';
end

end