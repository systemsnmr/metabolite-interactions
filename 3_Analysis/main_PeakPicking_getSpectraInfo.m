function h = main_PeakPicking_getSpectraInfo(spectra,fsamode,pathtofigs,ploton,localGlobalNoiseThreshold,peaks)

% Author Maren (based on Yaro's script)
%% Help 
% h = getPeaksAndSN2b_frompeaklist(spectrumPMP,spectrumM,fsamode,sncutoff,ploton,localGlobalNoiseThreshold,peaks,rf)
% 
% Takes 1D NMR spectrum, identifies peaks, calculates SN for every peak.
% Also calculates FSA parameter (or relaxation factor) - more quantitative representation of
% interaction characteristics. For this REFERENCE SPECTRA ARE REQUIRED
%
% Returns data structure with hits: point index of peak, ppm value, SN, intensity.
% (h.point, h.ppm, h.SN, h.intensity, h.FSA, h.FSA_freeM)
%
% NMR spectra are expected to be a text file exported from TopSpin:
% File > Save... (Ctrl+S) -> "save data of displayed region in a text file".
%
% Input arguments:
% spectra                    = structure with names and types of all nmr spectra
% fsamode                    = mode to calculate FSA (see main script)
% figureout                  = outputfolder for figures
% ploton                     = plot single peaks (1) or not (0)
% localGlobalNoiseThreshold  = when "noise" around the peak is high use it to determine final SN ratio.
% peaks                      = structure with info about peaks that can be found in the mix (peaks.left, peaks.right,peaks.Names)

%%%% reference spectra
% ref_1 - reference of M (free metabolite, may show FSA>1 due to Chem.Shift Perturbation of M signals.
% ref_2 - reference of PMP (Protein+Mix minus Protein. Should have no problem with CSP).

%% Set addiational parameters

% parameter to assess stability
stab_param = 0.05; % if S1/S2 > stab_param --> flag peak as not-stable


%% Error-check input arguments
if ~isnumeric(fsamode)
    error('fsamode has to be 0,1 or 2');
end
if ~isnumeric(ploton)
    error('ploton has to be 0 or 1');
end
if ~isnumeric(localGlobalNoiseThreshold)
    error('localGlobalNoiseThreshold has to be set');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Detect signals in the known peak windows - Maren: 03.03.2017

% get index of PMP and M spectrum
pmp_ind = strcmp([spectra.type],'T1rho_PMP');
m_ind = strcmp([spectra.type],'T1rho_M');

% get index of ref sepctra spectrum
ref1_ind = strcmp([spectra.type],'T1rho_200_M');
ref2_ind = strcmp([spectra.type],'T1rho_10_PMP');

% get index of spectra for calculation of relaxation factor
rf1_ind = strcmp([spectra.type],'T1rho_10_M');
rf2_ind = strcmp([spectra.type],'T1rho_200_M');
rf3_ind = strcmp([spectra.type],'T1rho_10_PMP');
rf4_ind = strcmp([spectra.type],'T1rho_200_PMP');

% get indices of S1 and S2 short spectra
s1_ind = strcmp([spectra.type],'T1rho_S1_M');
s2_ind = strcmp([spectra.type],'T1rho_S2_M');
   
% Spectrum import
for i = 1:length(spectra)
    try 
        [spectra(i).numPoints,spectra(i).leftPPM,spectra(i).rightPPM,spectra(i).intensity,spectra(i).ppmValues] = import_spectrum(spectra(i).fullfile);
    catch
        warning(strcat('Missing file:',spectra(i).fullfile));
    end
end

% Calc  default noise within 1 ppm (1*numberOfPoints/spectralWidth) of the start of the spectrum
% use PMP and M spectrum for noisecalc
globalNoise = noisecalc(spectra(pmp_ind),spectra(m_ind));


%%
numOfHits = length(peaks.right);

% Preallocate to store indices, intensities, noise and S/N
h = struct;
h(1).Index = []; % peak index (in the points of spectrum)
h(1).Intensity = [];
h(1).SNglobal = []; % S/N based on the global noise in the spectrum
h(1).SNfinal = []; % Final value of S/N selected for export
h(1).flagFSA = []; % to exclude values with neg. FSA (0: bad, 1: good);
h(1).flagSN = []; % to exclude values with neg. SN (0: bad, 1: good);
h(1).flagS1S2 = []; % to mark if metabolites are unstable in short term (0: no, 1: yes);
h(1).flagNeg = []; % to mark if metabolites are more stable in the presence of protein (0: no, 1: yes);

% Get the indices and intensities of the strongest signals in every hit.
% Calculate S/N, normalized S/N and pick the optimal one.
for i=1:numOfHits
    
    % Get index of max value, and the value itself:
    for j = 1:length(spectra)
        % find all peak regions in all spectra
        [spectra(j).peaks(i).LB, spectra(j).peaks(i).RB, spectra(j).peaks(i).Width] = find_peaks(spectra(j).ppmValues,peaks.left(i),peaks.right(i));
        [spectra(j).peaks(i).Index,spectra(j).peaks(i).ppm,spectra(j).peaks(i).Intensity] = ...
            find_max(spectra(j).intensity,spectra(j).peaks(i).LB,spectra(j).peaks(i).RB,spectra(j).ppmValues);
    end

    % save index and ppm of peak
    h(i).Index = spectra(pmp_ind).peaks(i).Index;
    h(i).ppm = spectra(pmp_ind).peaks(i).ppm;
    
    % calculate subtracted intensity:
    h(i).Intensity = spectra(pmp_ind).peaks(i).Intensity-spectra(m_ind).peaks(i).Intensity;
    
    % S/N in reference to general noise:
    h(i).SNglobal = h(i).Intensity/globalNoise;
    
    % caluclate local noise
    [h(i).LocalNoise, h(i).NoiseMedian] = calclocalnoise(spectra(pmp_ind),spectra(m_ind),i);

    % Calculate varios versions of noise 
    h(i).SN = h(i).Intensity / h(i).LocalNoise;
    % adding absolute, not to get negative medians
    h(i).CorrectedNoise = abs(h(i).LocalNoise - abs(h(i).NoiseMedian)); 
    % max - to avoid very large S/N when noise is very small
    h(i).SNM = h(i).Intensity / max([globalNoise,h(i).CorrectedNoise]);
    
    if h(i).LocalNoise/globalNoise < localGlobalNoiseThreshold
        h(i).SNfinal = h(i).SNglobal; % take global noise value
    else
        h(i).SNfinal = h(i).SNM; % take local noise, corrected with median
    end
    
    % flag if S/N < 0
    if h(i).SNfinal < 0
        h(i).flagSN = 0;
    else
        h(i).flagSN = 1;
    end
    
    % calculate different version of FSA: (PMP-M)/M
    if fsamode == 0
        h(i).FSA      = h(i).Intensity / spectra(m_ind).peaks(i).Intensity;
    elseif  fsamode == 1 % actual FSA
        h(i).FSA       = h(i).Intensity / spectra(ref2_ind).peaks(i).Intensity;
    elseif fsamode == 2 % relaxation factor
        h(i).RF_noprot   = (spectra(rf2_ind).peaks(i).Intensity/spectra(rf1_ind).peaks(i).Intensity);
        h(i).RF_plusprot = (spectra(rf4_ind).peaks(i).Intensity/spectra(rf3_ind).peaks(i).Intensity);
        h(i).FSA         = h(i).RF_noprot - h(i).RF_plusprot;
    end
    if h(i).FSA < 0
        h(i).flagFSA = 0;
    else
        h(i).flagFSA = 1;
    end

    % calculate FSA_free M
    h(i).FSA_freeM = h(i).Intensity / spectra(ref1_ind).peaks(i).Intensity;
    
    %% check metabolite stability
    h(i).flagS1S2 = 1;
    h(i).flagNeg = 1;
    
    % look at S1 and S2 spectra -> average over peak width to remove bias due to phase shift
    s1mean = mean(spectra(s1_ind).intensity(spectra(s1_ind).peaks(i).LB:spectra(s1_ind).peaks(i).RB));
    s2mean = mean(spectra(s2_ind).intensity(spectra(s2_ind).peaks(i).LB:spectra(s2_ind).peaks(i).RB));
    if abs(1-(s1mean/s2mean))>stab_param
         h(i).flagS1S2 = 0;
    end
    
    % check if metabolites are more stable in the presence of protein
    % results in NEGATIVE intensities: M_short - [PM_short - P_short] (= T1r10freeM_prot_mix.txt – T1r10PMP_prot_mix.txt)
    stab = spectra(rf1_ind).peaks(i).Intensity - spectra(rf3_ind).peaks(i).Intensity;
    if stab < -stab_param
        h(i).flagNeg  = 0;
    end

end

%% Plot peaks if selected    
if ploton
    plotpeaks(pathtofigs,spectra);
end

%% Assemble structure for export

lh = numel(h);

[ex(1:lh).point] = deal(h.Index);
[ex(1:lh).ppm] = deal(h.ppm);
[ex(1:lh).SN] = deal(h.SNfinal);
[ex(1:lh).intensity] = deal(h.Intensity);
[ex(1:lh).FSA] = deal(h.FSA);
[ex(1:lh).FSA_freeM] = deal(h.FSA_freeM);
[ex(1:lh).flagSN] = deal(h.flagSN);
[ex(1:lh).flagFSA] = deal(h.flagFSA);
[ex(1:lh).flagS1S2] = deal(h.flagS1S2);
[ex(1:lh).flagNeg] = deal(h.flagNeg);
clear h;
h = ex;
end % function


function [numPoints,leftPPM,rightPPM,intensity,ppmValues] = import_spectrum(spectrum)

dataFile = importdata(spectrum);

% Get number of points in the spectrum:
numPoints = strread(dataFile.textdata{6},'# SIZE = %d ( = number of points)');
% Get left & right boundaries of Spectral Region:
[leftPPM, rightPPM] = strread(dataFile.textdata{4},'# LEFT = %f ppm. RIGHT = %f ppm.');

ppmValues = linspace(leftPPM,rightPPM,numPoints)'; % Transposed!
intensity = dataFile.data; 

end

function globalNoise = noisecalc(PMP,M)

%check how different PMP and M spectrum are, error if too different
if abs(1-M.numPoints/PMP.numPoints) > 0.001
    error('Implement new routine here!');
end

% otherwise, calculate index from ppm= 14 -> 13
[~,indstart1] = min(abs([M.ppmValues]-13));
[~,indstart2] = min(abs([PMP.ppmValues]-13));
[~,indend1] = min(abs([M.ppmValues]-12));
[~,indend2] = min(abs([PMP.ppmValues]-12));

indstart = min(indstart1,indstart2);
indend = min(indend1,indend2);

% check if too close to start
if indstart < 5000
    error('Too close to start of spectrum!');
end

% get spectrum for PMP and M
mint = M.intensity(indstart:indend);
pmpint = PMP.intensity(indstart:indend);

subint = pmpint-mint;
globalNoise = std(subint);

end

function plotpeaks(pathtofigs,spectra)

   colormatrix = colormap(jet(4));
   groupn = {'PMP','M','T1r200freeM','T1r10PMP'};
   
   % get indices of spectra
   ind = [];
   ind(1) = find(strcmp([spectra.type],'T1rho_PMP'));
   ind(2) = find(strcmp([spectra.type],'T1rho_M'));
   ind(3) = find(strcmp([spectra.type],'T1rho_200_M'));
   ind(4) = find(strcmp([spectra.type],'T1rho_10_PMP'));
 
for i = 1:length(spectra(1).peaks)

    filename = char(strcat(pathtofigs,'\',num2str(i),'.png'));
    
    f = figure('visible','off');
    hold all;
    
    for j = 1:length(groupn)
        xvals = spectra(ind(j)).ppmValues(spectra(ind(j)).peaks(i).LB:spectra(ind(j)).peaks(i).RB);
        yvals = spectra(ind(j)).intensity(spectra(ind(j)).peaks(i).LB:spectra(ind(j)).peaks(i).RB);

        color = colormatrix(j,:);
        a(j) = plot(xvals,yvals,'Color',color,'LineWidth',2);    
        
        if j == 1
            miny = min(yvals);
            maxy = max(yvals);
            minx = min(xvals);
            maxx = max(xvals);
        end
    end

    pmpind =  spectra(ind(1)).peaks(i).ppm;
    mind = spectra(ind(2)).peaks(i).ppm;
    
    % plot lines
    line([pmpind pmpind] , [miny maxy],'Color', colormatrix(1,:)); 
    line([mind mind] , [miny maxy],'Color', colormatrix(2,:));

    legend([a], [groupn],'FontSize',10,'location', 'southeast');

    xlim([minx maxx]);
  
    saveas(f,filename);
    
    clearvars xvals yvals;
end

end
 

function [LocalNoise, NoiseMedian] = calclocalnoise(PMP,M,i)

% peakWidth                             ----
% sigWidth  (4*peakWidth)             --------
% noiseSpan (12*peakWdith)     ------------------------ 
% noiseWidth == sigWidth      --------

pmpinfo = struct;
minfo = struct;

% define region for local noise
pmpinfo.peakWidth = PMP.peaks(i).RB - PMP.peaks(i).LB;
minfo.peakWidth = M.peaks(i).RB - M.peaks(i).LB;

% check if PMP and M have very different widths
if abs(pmpinfo.peakWidth-minfo.peakWidth) > 10
   disp(strcat('PMP and M very different:',num2str(abs(pmpinfo.peakWidth-minfo.peakWidth)),'points')); 
end

pmpinfo.leftstart = PMP.peaks(i).LB - round(5.5*pmpinfo.peakWidth,0);
pmpinfo.leftend = PMP.peaks(i).LB - round(1.5*pmpinfo.peakWidth,0);
pmpinfo.rightstart = PMP.peaks(i).RB + round(1.5*pmpinfo.peakWidth,0);
pmpinfo.rightend = PMP.peaks(i).RB + round(5.5*pmpinfo.peakWidth,0);

minfo.leftstart = M.peaks(i).LB - round(5.5*minfo.peakWidth,0);
minfo.leftend = M.peaks(i).LB - round(1.5*minfo.peakWidth,0);
minfo.rightstart = M.peaks(i).RB + round(1.5*minfo.peakWidth,0);
minfo.rightend = M.peaks(i).RB + round(5.5*minfo.peakWidth,0);

% check that number of points is the same
if ~((minfo.leftstart-minfo.leftend) == (pmpinfo.leftstart-pmpinfo.leftend))
    difflen = abs((minfo.leftend-minfo.leftstart) - (pmpinfo.leftend-pmpinfo.leftstart));
    pmpinfo.addsub = round(difflen/2,0);
    if mod(difflen,2) > 0
        minfo.addsub = pmpinfo.addsub + 1;
    else
        minfo.addsub = pmpinfo.addsub;
    end
    if (minfo.leftend-minfo.leftstart) > (pmpinfo.leftend-pmpinfo.leftstart)
        minfo.leftstart = minfo.leftstart + minfo.addsub;
        pmpinfo.leftstart = pmpinfo.leftstart - pmpinfo.addsub;
    else
        minfo.leftstart = minfo.leftstart - minfo.addsub;
        pmpinfo.leftstart = pmpinfo.leftstart + pmpinfo.addsub;
    end
    
end
if ~((minfo.rightstart-minfo.rightend) == (pmpinfo.rightstart-pmpinfo.rightend))
    difflen = abs((minfo.rightend-minfo.rightstart) - (pmpinfo.rightend-pmpinfo.rightstart));
    pmpinfo.addsub = round(difflen/2,0);
    if mod(difflen,2) > 0
        minfo.addsub = pmpinfo.addsub + 1;
    else
        minfo.addsub = pmpinfo.addsub;
    end
    if (minfo.rightend-minfo.rightstart) > (pmpinfo.rightend-pmpinfo.rightstart)
        minfo.rightstart = minfo.rightstart + minfo.addsub;
        pmpinfo.rightstart = pmpinfo.rightstart - pmpinfo.addsub;
    else
        minfo.rightstart = minfo.rightstart - minfo.addsub;
        pmpinfo.rightstart = pmpinfo.rightstart + pmpinfo.addsub;
    end
    
end

pmpinfo.noiserange = [pmpinfo.leftstart:pmpinfo.leftend  pmpinfo.rightstart:pmpinfo.rightend];
minfo.noiserange = [minfo.leftstart:minfo.leftend  minfo.rightstart:minfo.rightend];


% compute difference spectra in noise range
intdiff = PMP.intensity(pmpinfo.noiserange)-M.intensity(minfo.noiserange);

LocalNoise = std(intdiff);    
NoiseMedian = median(intdiff);

end

function [LB, RB, Width] = find_peaks(ppmValues,pleft,pright)
    % find left indices
    [~,LB]= min(abs(ppmValues-pleft));
     % find right indices
    [~,RB]= min(abs(ppmValues-pright));
    Width = ppmValues(LB)-ppmValues(RB);
end

function [idx,ppm,maxint] = find_max(intensity,LB,RB,ppmValues)
    [maxint, index] = max(intensity(LB:RB));
    idx = LB+index-1;
    if nargin > 3
        ppm = ppmValues(idx);
    else
        ppm = [];
    end
end
