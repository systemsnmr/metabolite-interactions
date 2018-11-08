function main_PeakPicking(filePath,mix,protein,LGThreshold,fsamode,ploton,input,output)

% Author: Yaroslav Nikolaev. 2015-06-04.
% adapted by Maren Diether

%% Help / annotation
% call: getHitsAssignmentsAndExportResults_fromPeaklist(filePath,mix,protein,LGThreshold,fsamode,ploton,input,output)
%
% input
% filePath: path where PMP and M spectra can be found
% mix: name of mix to be analyzed
% protein: name of protein to be analyzed
% LGThreshold: Local Global Noise Threshold
% fsamode: determine fsa mode calculation (see main.m)
% ploton: 1: plot every single peak, 0: don't
% input: input files & folders
% output: output files & folders
%
% output
% hits are exported to results-folder (set in main.m)
%
%% Checking inputs
%==========================
if or(~ischar(mix),~ischar(protein))
    error('Mix and protein parameter should be a string');
end

if ~isempty(LGThreshold) && ~isnumeric(LGThreshold)
    error('Local Global Noise cutoff should be a numerical variable');
end

% set exptName
exptName = strcat(protein,'_',mix);
fprintf(1,'== Experiment: %s \n', exptName);

% make new structure for info about all spectra
spectra = struct;
for i=1:length(input.spectranames.type)
    spectra(i).type = input.spectranames.type(i);
    spectra(i).fullfile = generateFilename(protein,mix,spectra(i).type,input.spectranames,filePath);
end


%% Compound names, and PPM regions of compound peaks in the mix:
%=============================================================

% define current peaklist
peaklist = fullfile(input.peaklist,strcat(mix,'.csv'));

% check file
if exist(peaklist, 'file') ~= 2
    error('Peaklist file (%s) not found',peaklist);
end

% read peaklist
fileID = fopen(peaklist);
raw = textscan(fileID,'%*f,%f,%f,%s');
fclose(fileID);

peaks.left = raw{1}+0.33*raw{2};
peaks.right = raw{1}-0.33*raw{2};

% parse peakNames
% remove '"' and split by ','
for i = 1:length(raw{3})
    splitname = raw{3}{i,1}(2:end-1); % to remove the '"'
    ll = length(strsplit(splitname,','));
    peaks.names(i,1:ll) = strsplit(splitname,',');
end

% get unique list of compound Names
compoundNames = reshape(peaks.names,[],1);
compoundNames( all(cellfun(@isempty,compoundNames),2), : ) = []; % delete empty fields
compoundNames = unique(compoundNames);

clearvars raw splitname ll;

%% Load hits
%============================================================

% make folder for peak figures if selected
if ploton
    pathtofigs = char(fullfile(output.figures,protein,mix));
    pathsub =  char(fullfile(output.figures,protein));
    if ~(exist(pathsub,'dir'))
        mkdir(pathsub);
    end
    if ~(exist(pathtofigs,'dir'))
        mkdir(pathtofigs);
    end
else
    pathtofigs = '';
end

h = main_PeakPicking_getSpectraInfo(spectra,fsamode,pathtofigs,ploton,LGThreshold,peaks);

% output: h.ppm, h.SN, h.point, h.intensity, h.FSA, h.FSA_freeM,
% and flags for usability: h.flagSN, h.flagFSA, h.flagS1S2, h.flagNeg

%% Identify unique hits
%=============================================================
numOfHits = numel(h);
h(1).allMatches = {};
h(1).assign = {};
h(1).unique = [];
h(1).final = [];

for i=1:numOfHits
    % Keep all matches as cell array - easier to visualize ambigs later:
    currnames = peaks.names(i,:);
    currnames(:,find(all(cellfun(@isempty,currnames),1))) = []; %delete whole empty column
    h(i).allMatches = currnames;
    if size(currnames,2)>1 % if more than one compound match found
        h(i).assign = strjoin(currnames);
        h(i).unique = 0; % Set hit as non-unique.
    else
        h(i).assign = currnames{1};
        h(i).unique = 1; % Set hit as unique.
    end
    h(i).final = 0;
end

%% Export results
%=============================================================
clear hits;

hits = dataset({[1:numOfHits]','hit'},...
    {[h.point]','point'},...
    {[h.ppm]','PPM'},...
    {[h.SN]','SN'},...
    {[h.FSA]','FSA'},...
    {[h.FSA_freeM]','FSA_freeM'},...
    {{h.allMatches}','AllMatches'},...
    {[h.intensity]','Intensity'},...
    {[h.unique]','Unique'},...
    {{h.assign}','Assign'},...
    {[h.final]','Final'},...
    {[h.flagFSA]','FlagFSA'},...
    {[h.flagSN]','FlagSN'},...
    {[h.flagS1S2]','FlagS1S2'},...
    {[h.flagNeg]','FlagNeg'});

save(strcat(output.hits, 'T1r_',exptName, '.mat'), 'hits');


end % fullAnalysis

function fname = generateFilename(protein,mix,type,namesinfo,fpath)

% figure out index
idx = find(strcmp(namesinfo.type,type));

fname = strcat(namesinfo.prefix(idx),'_',protein,'_',mix);
if ~isempty(namesinfo.suffix{idx})
    fname = strcat(fname,'_',namesinfo.suffix(idx),'.txt');
else
    fname = strcat(fname,'.txt');
end

if nargin > 4
    fname = fullfile(fpath,fname);
end

fname = char(fname);

end
