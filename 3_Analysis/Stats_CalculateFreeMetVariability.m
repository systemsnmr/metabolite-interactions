function Stats_CalculateFreeMetVariability

% function to calculate the variability of signal in free metabolite mixes
%
% input
% the script expects the mix-spectra to be in 'pathtoM'
% the peaklists are expected to be in the 'peakfolder'
%
% output
% results figures are saved in 'resfolder'

%% set params
mixes = {'ccm1','ccm2','ccm3','ccm4'};

% path to all free M spectra
pathtoM = '.\input_spectra\';

% path to peaklists
peakfolder = '.\input_peaklists\';

% folder for results
resfolder = '.\results_metabolitestats\';

%%
% get all free M spectra
specs = get_freeMspectra(pathtoM,mixes);

%% perform analysis for all peaklists
for i = 1:length(mixes)
   
    % name of current peaklist
    peaklistname = char(strcat(peakfolder,mixes(i),'.csv'));
    
    % calculate intensities of each peak in each spectrum
    [h,pnames] = get_intensities(specs.(mixes{i}),peaklistname);
    
    % plot statitics for this mix
    plotname = char(strcat(resfolder,mixes(i),'.png'));
    plot_results(plotname,h,pnames);
   
    % clear results structure
    clear h;
    
end

end 


%% sub-functions
function [h,pnames] = get_intensities(spectrumM,peakfile)
% calculate intensities of each peak in each spectrum

% peakfile import
peaks = import_peaklist(peakfile);
pnames = peaks.names;

% Spectrum import
M = struct;
for i = 1:length(spectrumM) 
    [M(i).numPoints,M(i).leftPPM,M(i).rightPPM,M(i).intensity,M(i).ppmValues] = import_spectrum(spectrumM{i});
end
    
h = struct;

for i = 1:length(peaks.right) 

    for j = 1:length(spectrumM)
        % find left indices
        [~,h(j).LBM(i)]= min(abs(M(j).ppmValues-peaks.left(i)));
         % find right indices
        [~,h(j).RBM(i)]= min(abs(M(j).ppmValues-peaks.right(i)));
        h(j).WidthM(i) = M(j).ppmValues(h(j).LBM(i))-M(j).ppmValues(h(j).RBM(i));
        
        % get max values
        h(j).maxValue(i) = max(M(j).intensity(h(j).LBM(i):h(j).RBM(i)));
    end
    
end


end 

function [numPoints,leftPPM,rightPPM,intensity,ppmValues] = import_spectrum(spectrum)

% import NMR spectra from txt files

dataFile = importdata(spectrum);

% Get number of points in the spectrum:
numPoints = strread(dataFile.textdata{6},'# SIZE = %d ( = number of points)');
% Get left & right boundaries of Spectral Region:
[leftPPM, rightPPM] = strread(dataFile.textdata{4},'# LEFT = %f ppm. RIGHT = %f ppm.');

ppmValues = linspace(leftPPM,rightPPM,numPoints)'; % Transposed!
intensity = dataFile.data; 

end


function peaks = import_peaklist(peakfile)

% import the peaklists from .csv files

% read file
if exist(peakfile, 'file') ~= 2
    error('Peaklist file (%s) not found',peakfile);
end

% read peaklist
fileID = fopen(peakfile);
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

end


function plot_results(plotname,h,pnames)
% plot peak intensities for every mix

% find peaks that have unique hits
idx = [];
for i = 1:size(pnames,1)
    if isempty(pnames{i,2})
        if ~strcmp(pnames{i,1},'buffer')
            idx(end+1) = i;
        end
    end
end


% assemble x-axis labels from peakindex & compound name
pnames = pnames(idx);
for i = 1:length(pnames)
   pnames(i) = strcat(num2str(idx(i)),'_',pnames(i));
end

% assemble matrix for boxplot
plotmat = [];
for j = 1:length(h)
   plotmat(j,:) = h(j).maxValue(idx); 
end

% plot
f = figure('visible','off');
set(f, 'Position', [0 100 1650 1000]); % [x y width height]

% divide in 2 subplots
% 1: boxplot of intensities
plot_handle = subplot(2,1,1);
hold all;
boxplot(plotmat,'PlotStyle','compact','Labels',pnames);
xlabel('compounds');
ylabel('intensity');
title('Intensity of unique peaks');
hold off;

% 2: plot of variability
plot_handle = subplot(2,1,2);

hold all;
varmat = std(plotmat)./mean(plotmat);
varmat(2,:) = NaN;
boxplot(varmat,'PlotStyle','compact','Labels',pnames);
ylim([0 0.7]); % fix so that it is the same everywhere
xlabel('compounds');
ylabel('Std(peaks)/mean(peaks)');
title('Variability');
hold off;

saveas(f,plotname);

end

function specs = get_freeMspectra(indir,mixes)

acontent = dir(indir);
allcontent = {acontent(3:end).name};

specs = struct;
for m = 1:length(mixes)
    specs.(mixes{m}) = {};
end

% read content of folder and extract pnames and mixnames from there
for i = 1:length(allcontent)
    
    % get path to txts
    currpath = char(fullfile(indir,allcontent(i)));
    content = dir(currpath);
    
    % split filenames
    allfnames = {content.name};
    allfnames = allfnames(strncmp('T1r200freeM_',allfnames,10));
    
    for j = 1:length(allfnames)
        for m = 1:length(mixes)
            if find(regexpi(allfnames{j},char(strcat('T1r200freeM_\w{3,4}_',mixes(m),'.txt'))))
                specs.(mixes{m})(end+1) = fullfile(currpath,allfnames(j));
            end
        end
    end
end
    
end