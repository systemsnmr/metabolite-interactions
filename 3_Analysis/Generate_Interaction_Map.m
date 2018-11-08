function [stats,td] = Generate_Interaction_Map(zscores,lower_cutoff,sncutoff,database,scoringopt,special)

%%
% Input
% zscores:      0: do not convert map to zscores
%               1: convert map to zscores
% lower cutoff: any value between 0 and 1:  set all FSI/RF lower than this to 0 (only when z-scores = 0)
% sncutoff:     any value > 0:              all peaks with S/N lower than that will be ignored
% database:     'e': only ecocyc
%               'eb': ecocyc & brenda database
% scoringopt:   1: filter peaks based on lower cutoff and s/n cutoff separately
%               2: filter peaks based on the product of lower cutoff and s/n cutoff 
% special:      0: normal run
%               1: fast mode - if hit structure is provided here, this is used (useful for fast parameter scanning)
%
% Output
% Various figures, statistics and tables
% In the output folder
%
% Example call
% Generate_Interaction_Map(0,0.1805,2,'e',1)
%
% ============================================================
%% Beginning of script
%

% check input
if nargin < 6
    special = [];
end
if nargin < 5
    % use default parameters
    database = 'e';
    sncutoff = 2;
    lower_cutoff = 0.1805;
    zscores = 0;
    scoringopt = 1;
end
if zscores == 1
    lower_cutoff = 0;
end

%% more parameters and variables
medianopt = 1; % compute median peak score of all peaks per compound that survive filtering

% folder for output
stats = struct;
results_folder = './results_tmp';
if ~exist(results_folder, 'dir')
    mkdir(results_folder); % folder for output
    disp('Created folder for figure output');
end

% input data
input = struct;
input.data = './results_matfiles/';
input.files = dir(strcat(input.data,'*.mat') );

% input on compound similarities
cpdsim = load('input_other\Simcomp2_compoundsimilarities.mat');
input.cpdsim = cpdsim.output;
clear cpdsim;

% load database info on known interactions
if find(strcmp(database,'e'))
    load('input_other\knownInteractions.mat');
elseif find(strcmp(database,'eb'))
    load('input_other\knownInteractions_withBrenda.mat');
else
   error('database not specified'); 
end

% plotting
fig = struct;
fig.font = 11;                      % font size in plot
fig.size = 30;                      % size of plot
absolute_colormap_threshold = 0;    %0: no threshold
                                    %any other value: set this as absolute threshold
% output
results_data_folder = char(strcat(results_folder,'/add_data'));
if ~exist(results_data_folder, 'dir')
    mkdir(results_data_folder); % folder for output
    disp('Created folder for data output');
end

% make structure for saving info on recovered interactions
stats.detected = 0;
stats.notdetected = 0;
stats.rdetected = 0;
stats.rnotdetected = 0;
stats.sdetected = 0;
stats.snotdetected = 0;
stats.overview = struct;
stats.protSize = [];

%% run analysis
if isempty(special)
    disp('Identifying hits...');
    [td.hUnfiltered, td.hUnique, td.mNames, td.pNames, addinfo] = regenerate_full_dataset(input.data,input.files,results_data_folder);
        
    disp('Generate Map...');
    [stats,mapFull] = generate_map(td.hUnfiltered,td.hUnique, td.pNames, td.mNames, zscores,lower_cutoff,sncutoff,knownInteractions,stats,medianopt,results_folder,results_data_folder,scoringopt);

    disp('Generate Figures...');
    Generate_Interaction_Map_MakeFigures(absolute_colormap_threshold,fig.font,fig.size,zscores,mapFull,td.pNames,td.mNames,knownInteractions,stats,addinfo.mix,results_folder,results_data_folder);
    
    disp('Generate additional figures...');
    Generate_Interaction_Map_MakeMetaboliteOverview(mapFull,td.pNames,td.mNames,knownInteractions,results_folder);
    Generate_Interaction_Map_ComputeCompoundSimilarities(mapFull,input.cpdsim,'local',td.pNames,td.mNames,knownInteractions,results_folder);
    Generate_Interaction_Map_ComputeCompoundSimilarities(mapFull,input.cpdsim,'global',td.pNames,td.mNames,knownInteractions,results_folder);
    Generate_Interaction_Map_MakeOutliersFigures(fig.font,fig.size,addinfo, td.pNames,td.mNames,stats,results_folder);
    
else % fast mode for parameter scanning
    [stats,~] = generate_map(special.hUnfiltered, special.hUnique, special.pNames, special.mNames, zscores,lower_cutoff, sncutoff, knownInteractions,stats,medianopt,'','',scoringopt);
    td = [];
end
end

%============================================================
function [hUnfiltered, hits, mNames, pNames, addinfo]= regenerate_full_dataset(data_folder,files, results_folder)
%============================================================

%% load all info on peaks and combine in one big 'hit' structure

n_files = size(files,1);
h = cell(n_files,1);

for i=1:n_files
    tmp = load( strcat(data_folder,files(i).name) );
    h{i} = tmp.hits;
    clear tmp;
    
    n_hits = size(h{i},1);
    
    nameParse = textscan(files(i).name(1:end-4),'%s','delimiter','_'); % splits by underscore delimeter   
    [expt, prot, mix] = deal(nameParse{1}(1), nameParse{1}(2), nameParse{1}(3));

    % change name if protein is RpiA or TpiA (typo)
    if strcmp(prot,'RpiA')
        prot = {'rpia'};
    end
    if strcmp(prot,'TktA')
        prot = {'tkta'};
    end
    
    h{i}.protein = repmat(prot, n_hits, 1);
    h{i}.experiment = repmat(expt, n_hits, 1);
    h{i}.mix = repmat(mix, n_hits, 1);
    
    %%% Merge datasets
    if i==1
        hits = h{1};
    else
        hits = cat(1,hits,h{i});
    end
end

% convert to structure
hits = dataset2struct(hits);

%% filtering
% keep only unique hits
hits = hits([hits.Unique]==1);

% remove buffer
hits = hits(~strcmp({hits.Assign},'buffer'));

% remove proteins not included in this study
rmnames = {'BSA','lrp','cra','fis'};
for i = 1:length(rmnames)
    hits = hits(~strcmp({hits.protein},rmnames{i}));
end

% define new column as peak-id
for i = 1:length(hits)
    hits(i).peakid = strcat(hits(i).Assign,'-',num2str(hits(i).hit));
end

%%
addinfo = struct;
addinfo.hNeg =   hits([hits.FlagNeg]==0);
addinfo.hS1S2 = hits([hits.FlagS1S2]==0);

mix = unique({hits.mix});
for i = 1:length(mix)
    currh = hits(strcmp({hits.mix},mix(i)));
    peakNames(i).Names = unique({currh.peakid});
end
addinfo.peakNames = peakNames;

% plot statistics on peaks (before filtering)
plot_peak_stats(hits,results_folder);

% save hit structure
hUnfiltered =  hits;

%% filtering part II

% set all FSI/RF below lower cutoff to 0 and remove non-optimal peaks
idx = [];
for i = 1:length(hits)
    if sum([hits(i).FlagNeg hits(i).FlagS1S2 hits(i).FlagFSA hits(i).FlagSN]) < 4
        idx(end+1) = i;
    end
end
hits(idx) = [];

% find compound names
mNames = unique({hits([hits.Unique]==1).Assign});

% find mix association for every metabolite
metmix = [];
for j = 1:length(mNames)
    mmh = hits(strcmp({hits.Assign},mNames(j)));
    metmix(j) = str2num(mmh(1).mix(4));
end
addinfo.mix = metmix;
    
% find protein names
pNames = unique({hits.protein});

% save the merged hit structure
if ~exist(fullfile(data_folder,'merged_data'), 'dir')
  mkdir(fullfile(data_folder,'merged_data'));
end

save(fullfile(data_folder, 'merged_data', 'all_hits.mat'), 'hits');

end


%============================================================
function [stats,mapFull] = generate_map(hits,hUnique, pNames, mNames, zscores,lower_cutoff,sncutoff,knownInteractions,stats,medianopt,results_folder,results_data_folder,scoringopt)
%============================================================

hU = hUnique;
idx = [];

% do the last part of filtering:
if scoringopt == 1
    for i = 1:length(hU)
        if hU(i).FSA < lower_cutoff
            idx(end+1) = i;
        elseif hU(i).SN < sncutoff
            idx(end+1) = i;
        end
    end
elseif scoringopt == 2
    % determine peak ranks
    allFSA = sort(unique([hU.FSA]),'ascend'); 
    allSN = sort(unique([hU.SN]),'ascend'); 
    maxrankFSA = length(allFSA);
    maxrankSN = length(allSN);
    maxrank = geomean([maxrankFSA maxrankSN]); % for scaling
    
    for i = 1:length(hU)
        rankFSA = 1 + maxrankFSA - mean(find(allFSA == hU(i).FSA));
        rankSN =  1 + maxrankSN - mean(find(allSN == hU(i).SN));
        
        currrank = geomean([rankFSA rankSN])/maxrank ;
        
        if currrank > lower_cutoff
            idx(end+1) = i;
        end
    end            
else
    error('wrong input');
end

hU(idx) = []; % delete peaks

% mixnames
mix = unique({hU.mix});

%% Generate hitMap & save as xlsx
if isempty(results_folder)
    filename = '';
else
    filename = fullfile(results_folder,'FinalResults.xlsx');
end
if medianopt
    % Generate hitMap from medians 
    [mapFull,snmap] = ComputeHitMap(pNames,'protein',mNames,'Assign',hU,zscores,filename,1);
else
    [mapFull,snmap] = ComputeHitMap(pNames,'protein',mNames,'Assign',hU,zscores,filename,0);
end

%% Generate hitMap for single mixes and save as xlsx
if ~isempty(results_folder)
    for i = 1:length(mix)
       currh = hU(strcmp({hU.mix},mix(i))); 
       currpeaks = unique({currh.peakid});
       filename = char(strcat(results_data_folder,'\',mix(i),'_Results.xlsx'));
       ComputeHitMap(pNames,'protein',currpeaks,'peakid',currh,0,filename,0);
       
       % compare density of peaks with original density
       filename = char(strcat(results_data_folder,'\',mix(i),'_DensityPlot.png'));
       CompareDensity(hits, currh, mix(i), length(pNames), filename);
    end
end

%% compute statitics on full FSI matrix

rmapFull = reshape(mapFull,1,[]);
numpos = length(find(rmapFull>0));

% stats on all hits
stats.numofvalues = length(rmapFull);
stats.valabovezero = numpos;
stats.valabovezeroperc = numpos/length(rmapFull);


% stats on known interactions
ctr = 1;
for i=1:length(pNames)
    pidx = find(strcmpi(knownInteractions.prot,pNames(i)));
    for t = 1:knownInteractions.oligomeric(pidx)
        stats.protSize(i,t) = knownInteractions.proteinSize(pidx)/1000;
    end
    
    for j = 1:length(mNames)
        % check if interaction is known and count
        midx = find(strcmpi(knownInteractions.met(pidx,:),mNames(j)));
        if ~isempty(midx)
            sflag = 0;
            rflag = 0;
            stats.overview(ctr).protein =  pNames(i);
            stats.overview(ctr).metabolites =  mNames(j);
            stats.overview(ctr).type =  knownInteractions.type(pidx,midx);
            stats.overview(ctr).protSize = knownInteractions.proteinSize(pidx)*knownInteractions.oligomeric(pidx);
            if strfind(stats.overview(ctr).type{1},'S')
                sflag = 1;
            end
            if strfind(stats.overview(ctr).type{1},'R')
                rflag = 1;
            end
            if ~isnan(knownInteractions.KMKI(pidx,midx))
                stats.overview(ctr).KMKI =  knownInteractions.KMKI(pidx,midx);
            else
                stats.overview(ctr).KMKI =  0;
            end
            stats.overview(ctr).FSI = mapFull(i,j);
            stats.overview(ctr).SN = snmap(i,j);
            if mapFull(i,j) > 0
                stats.overview(ctr).detected = 1;
                stats.detected = stats.detected+1;
                if sflag
                    stats.sdetected = stats.sdetected+1;
                end
                if rflag
                    stats.rdetected = stats.rdetected+1;
                end
            else
                stats.overview(ctr).detected = 0;
                stats.notdetected = stats.notdetected+1;
                if sflag
                    stats.snotdetected = stats.snotdetected+1;
                end
                if rflag
                    stats.rnotdetected = stats.rnotdetected+1;
                end
            end
            ctr = ctr +1;
        end
    end
end

stats.totalknown = stats.notdetected+stats.detected;
stats.detectedperc = stats.detected/stats.totalknown;


end

%%
function [mapFull,snmap] = ComputeHitMap(pNames,pfield,mNames,mfield,h,zscores,filename,medianopt)

mapFull = zeros(length(pNames),length(mNames)); % for FSIs
idxmap = cell(length(pNames),length(mNames)); % for indices of peaks
ppmmap = cell(length(pNames),length(mNames)); % for ppm
snmap = zeros(length(pNames),length(mNames)); % for S/N of peaks

if medianopt
    % compute median per compound
    for pIdx = 1:length(pNames) 
       for mIdx = 1:length(mNames)
           % find all matches
           allpids = strcmp({h.(pfield)},pNames(pIdx));
           allmids = strcmp({h.(mfield)},mNames(mIdx));
           allids = find((allpids+allmids)==2);
           if ~isempty(allids)
              mapFull(pIdx,mIdx) = median([h(allids).FSA]);
              snmap(pIdx,mIdx) =  median([h(allids).SN]);
              idxmap{pIdx,mIdx} = mat2str([h(allids).hit]); % also save idx of hit that was used for map
              ppmmap{pIdx,mIdx} = mat2str(round([h(allids).PPM],4));
           end
       end
    end    
else
    for i=1:length(h)
        pIdx = strcmp(pNames,h(i).(pfield));
        mIdx = strcmp(mNames,h(i).(mfield));
        
        mapFull(pIdx,mIdx) = h(i).FSA; % was SN here before!
        snmap(pIdx,mIdx) =  h(i).SN;
        if (sum(pIdx)+sum(mIdx)) > 1
            idxmap{pIdx,mIdx} = mat2str([h(i).hit]); % also save idx of hit that was used for map
            ppmmap{pIdx,mIdx} = mat2str([h(i).PPM]);
        end
    end    
end

%%% convert hitmap to z-scores if selected
if zscores == 1
    mapFull = zscore(mapFull,[],1);
end

%%% save indexes & FSAs to xlsx
if ~isempty(filename)
    sheet1 = 'Hit_indexes';
    sheet2 = 'FSI_matrix';
    sheet3 = 'SN_matrix';
    sheet4 = 'ppm';
    idxcell = {};
    idxcell(1,2:length(mNames)+1) = mNames;
    idxcell(2:length(pNames)+1,1) = pNames;
    idxcell(2:length(pNames)+1,2:length(mNames)+1) = idxmap;
    xlswrite(filename,idxcell,sheet1);
    idxcell(2:length(pNames)+1,2:length(mNames)+1) = num2cell(mapFull);
    xlswrite(filename,idxcell,sheet2);
    idxcell(2:length(pNames)+1,2:length(mNames)+1) = num2cell(snmap);
    xlswrite(filename,idxcell,sheet3);
    idxcell(2:length(pNames)+1,2:length(mNames)+1) = ppmmap;
    xlswrite(filename,idxcell,sheet4);
end

end

function plot_peak_stats(allpeaks,results_folder)

pout = struct;
poutmat = [];
mout = struct;

% aggregate per protein
pNames = unique({allpeaks.protein});
for i = 1:length(pNames)
   pout(i).name = pNames(i);
   cpeaks = allpeaks(strcmp({allpeaks.protein},pNames{i}));
   len = length(cpeaks);
   pout(i).FSA = (len-sum([cpeaks.FlagFSA]))/len;
   pout(i).Neg = (len-sum([cpeaks.FlagNeg]))/len;
   pout(i).S1S2 = (len-sum([cpeaks.FlagS1S2]))/len;
   poutmat(i,1) = pout(i).FSA;
   poutmat(i,2) = pout(i).Neg;
   poutmat(i,3) = pout(i).S1S2;
end

% aggregate per compound
mNames = unique({allpeaks.Assign});
for i = 1:length(mNames)
   mout(i).name = mNames(i);
   cpeaks = allpeaks(strcmp({allpeaks.Assign},mNames{i}));
   len = length(cpeaks);
   mout(i).FSA = (len-sum([cpeaks.FlagFSA]))/len;
   mout(i).Neg = (len-sum([cpeaks.FlagNeg]))/len;
   mout(i).S1S2 = (len-sum([cpeaks.FlagS1S2]))/len;
   moutmat(i,1) = mout(i).FSA;
   moutmat(i,2) = mout(i).Neg;
   moutmat(i,3) = mout(i).S1S2;
end

% plot stacked barplot
f = figure('visible','off');
set(f, 'Position', [0 100 1650 1000]); % [x y width height]
hold all;
bar(categorical(pNames),poutmat);
xtickangle(45);
ylim([0 1]);
ylabel('errors');
legend('FSA below 0','Stabilized by protein','degradation/conversion');
set(gca,'FontSize',16);
hold off;
fname = fullfile(results_folder,'peaks-proteins.png');
saveas(f,fname);

f = figure('visible','off');
set(f, 'Position', [0 100 1650 1000]); % [x y width height]
hold all;
bar(categorical(mNames),moutmat);
xtickangle(45);
ylim([0 1]);
ylabel('errors');
legend('FSA below 0','Stabilized by protein','degradation/conversion');
set(gca,'FontSize',16);
hold off;
fname = fullfile(results_folder,'peaks-metabolites.png');
saveas(f,fname);

% save
%pout= struct2table(pout);
%saveas(pout,'pout.xlsx');
%mout= struct2table(mout);
%saveas(mout,'mout.xlsx');
end

function CompareDensity(hits, currh, mix, lengthp, filename)
  
       allcurrh = hits(strcmp({hits.mix},mix)); 
       allh = unique([allcurrh.hit]);
       
       allpeaks = struct;
       for i = 1:length(allh)
           allpeaks(i).idx = allh(i);
           allpeaks(i).ppm = mean([allcurrh([allcurrh.hit]==allh(i)).PPM]); 
           % now, checked how often this ppm has been picked
           allpeaks(i).picked = length(currh([currh.hit]==allh(i))); 
           allpeaks(i).pickedpercent =  allpeaks(i).picked/lengthp;
           allpeaks(i).percent =  1;
       end
       
       f = figure('visible','off');
       hold all;
       set(f, 'Color', [1,1,1], 'Position', [100 100 1000 200]); % [x y width height]
       stem([allpeaks.ppm],[allpeaks.percent],'LineWidth',1.5,'Marker','none','Color',[0.8 0.8 0.8]);
       stem([allpeaks.ppm],[allpeaks.pickedpercent],'LineWidth',1.5,'Marker','none','Color',[0.22 0.62 0]);
       xlim([0 10]); 
       set(gca, 'XDir','reverse');
       ylabel('picked [%]');
       xlabel('ppm');
       title(['Picked Peaks in ' mix{1}]);
       hold off;
       saveas(f,filename);
       
end
       