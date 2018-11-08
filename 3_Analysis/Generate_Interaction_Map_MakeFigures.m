function Generate_Interaction_Map_MakeFigures(absolute_colormap_threshold,fSize,mapCellSize,zscores,mapFull,pNames,mNames,knownInteractions,stats,mixes,results_folder,results_data_folder)

%%% Display HitMap
fnamebase = fullfile(results_folder,'all_hits');
DisplayHitMap(absolute_colormap_threshold,fSize,mapCellSize,zscores,mapFull,pNames,mNames,knownInteractions,stats,fnamebase);

%%% Display HitMap sorted by mixes
fnamebase = fullfile(results_folder,'all_hits_sorted');
DisplayHitMap(absolute_colormap_threshold,fSize,mapCellSize,zscores,mapFull,pNames,mNames,knownInteractions,stats,fnamebase,mixes);

%% plot summary on interactions per protein / metabolite
binmapfull = mapFull;
binmapfull(mapFull>0) = 1;

results_folder2 = char(strcat(results_folder,'/stats_known_interactions/'));
if ~exist(results_folder2, 'dir')
    mkdir(results_folder2); % folder for output
    disp('Created folder for interaction output');
end

results_folder3 = char(strcat(results_folder,'/overview_detected_interactions/'));
if ~exist(results_folder3, 'dir')
    mkdir(results_folder3); % folder for output
    disp('Created folder for overview output');
end


% first: metabolites
plot_overview_interactions(binmapfull,1,mNames,results_folder3,'metabolite');

% second: proteins
plot_overview_interactions(binmapfull,2,pNames,results_folder3,'protein');


plot_percentual_interactions(binmapfull,mNames, pNames,results_folder2,knownInteractions);

%% plot stats

cmat = [0,0,1;1,0,0];

%%% save protein size versus # of hits
f = figure('visible','off');
hold all;

set(gca,'FontSize', 14,'LineWidth',2);

sizes = sum(stats.protSize,2);
ints = [];
for i = 1:size(mapFull,1)
    ints(i) = sum(mapFull(i,:)>0);
end

stem(sizes,ints,'filled');
%scatter(sizes,ints,[],cmat(1,:));

[c1,R1] = fitlnmodel(sizes,ints');
cnames = strcat('R2=',num2str(R1));

b(1) = plot(sizes,c1,'Color',cmat(2,:));

set(gca,'tickdir','out');
xlabel('protein size [kDa]');
ylabel('total #');
legend([b],[cnames]);
title('intensity vs. size');

%ylim([0 ylimm]);
hold off;
fname = fullfile(results_data_folder,'sizevs#.png');saveas(f,fname);



%%% 2 histogram of FSIs sorted by substrates and regulators
% S
IndexS = strfind([stats.overview.type], 'S');
idxS = find(not(cellfun('isempty', IndexS)));
% R
IndexR = strfind([stats.overview.type], 'R');
idxR = find(not(cellfun('isempty', IndexR)));

names = {'Substrates/Products', 'Regulators'};

f = figure('visible','off');
hold all;
a(1) = histogram(sort([stats.overview(idxS).FSI]),'FaceColor',cmat(1,:));
a(2) = histogram(sort([stats.overview(idxR).FSI]),'FaceColor',cmat(2,:));
if zscores == 1
    a(1).BinWidth = 0.5;
    a(2).BinWidth = 0.5;
else
    a(1).BinWidth = 0.025;
    a(2).BinWidth = 0.025;
end
legend([a],[names]);
xlabel('FSI');
ylabel('counts');
title('Histogram of FSI values of known interactions');
hold off;
fname = fullfile(results_data_folder,'Histogram_FSI.png');
saveas(f,fname);

%%% figure 3: FSI vs S/N
f = figure('visible','off');
hold all;
idx = find([stats.overview(idxS).FSI]>0);
a1(:,1) = [stats.overview(idxS(idx)).FSI];
b1(:,1) = [stats.overview(idxS(idx)).SN];

idx = find([stats.overview(idxR).FSI]>0);
a2(:,1) = [stats.overview(idxR(idx)).FSI];
b2(:,1) = [stats.overview(idxR(idx)).SN];

b(1) = scatter(a1,b1,[],cmat(1,:));
b(2) = scatter(a2,b2,[],cmat(2,:),'+');

% fit linear regression
[c1,R1] = fitlnmodel(a1,b1);
[c2,R2] =fitlnmodel(a2,b2);

cnames = [names, strcat('R2=',num2str(R1)), strcat('R2=',num2str(R2))];

b(3) = plot(a1,c1,'Color',cmat(1,:));
b(4) = plot(a2,c2,'Color',cmat(2,:));

xlabel('FSI');
ylabel('S/N');
title('FSI vs. S/N of known interactions');

if zscores == 1
    set(gca,'YAxisLocation','origin');
    xlimm = max(abs([a1;a2]));
    xlim([-xlimm xlimm]);
    legend([b],[cnames],'Location','northwest');
else
    legend([b],[cnames]);
end
ylimm = max(abs([b1;b2]));
ylim([0 ylimm]);
hold off;
fname = fullfile(results_data_folder,'SNvsFSI.png');
saveas(f,fname);

%%% figure 4: FSI vs KD of known interactions
f = figure('visible','off');

% find hits that have a KD
idxK = find([stats.overview.KMKI]>0);
idxRK = intersect(idxK,idxR);
idxSK = intersect(idxK,idxS);

hold all;
a1 = [];
a2 = [];
b1 = [];
b2 = [];
a1(:,1) = [stats.overview(idxSK).FSI];
a2(:,1) = [stats.overview(idxRK).FSI];
b1(:,1) = [stats.overview(idxSK).KMKI];
b2(:,1) = [stats.overview(idxRK).KMKI];

b(1) = scatter(b1,a1,[],cmat(1,:));
b(2) = scatter(b2,a2,[],cmat(2,:),'+');

cnames = names;

xlabel('KD');
ylabel('FSI');
title('FSI vs. KD of known interactions');

if zscores == 1
    set(gca,'XAxisLocation','origin');
    ylimm = max(abs([a1;a2]));
    %ylim([-ylimm ylimm]);
    legend([b],[cnames],'Location','northwest');
else
    legend([b],[cnames]);
end
set(gca, 'XScale', 'log');

hold off;
fname = fullfile(results_data_folder,'KDvsFSI.png');
saveas(f,fname);


%%% figure 5: FSI vs protein size
f = figure('visible','off');
hold all;
a1 = [];
a2 = [];
b1 = [];
b2 = [];

idx = find([stats.overview(idxS).FSI]>0);
a1(:,1) = [stats.overview(idxS(idx)).FSI];
b1(:,1) = [stats.overview(idxS(idx)).protSize];

idx = find([stats.overview(idxR).FSI]>0);
a2(:,1) = [stats.overview(idxR(idx)).FSI];
b2(:,1) = [stats.overview(idxR(idx)).protSize];

b(1) = scatter(a1,b1,[],cmat(1,:));
b(2) = scatter(a2,b2,[],cmat(2,:),'+');

% fit linear regression
[c1,R1] = fitlnmodel(a1,b1);
[c2,R2] =fitlnmodel(a2,b2);

cnames = [names, strcat('R2=',num2str(R1)), strcat('R2=',num2str(R2))];

b(3) = plot(a1,c1,'Color',cmat(1,:));
b(4) = plot(a2,c2,'Color',cmat(2,:));

xlabel('FSI');
ylabel('Protein [kDa]');
title('FSI vs. Protein size of known interactions');

if zscores == 1
    set(gca,'YAxisLocation','origin');
    xlimm = max(abs([a1;a2]));
    xlim([-xlimm xlimm]);
    legend([b],[cnames],'Location','northwest');
else
    legend([b],[cnames]);
end
hold off;
fname = fullfile(results_data_folder,'kDavsFSI.png');
saveas(f,fname);


end

function compute_colormap(zscores,mapFull,threshold)
startColor = [0.84 1 0.76]; % should not be completely white, otherwise hits with low S/N are almost not visible
endColor = [0.22 0.62 0];
numOfGradations = 1000;

% find lowest and highest value
minval = abs(min(min(mapFull)));
maxval = abs(max(max(mapFull)));

% compute gradient in pos and neg space
neggrad = round(numOfGradations*minval/(minval+maxval),0);
posgrad = round(numOfGradations*maxval/(minval+maxval),0);

if zscores == 1
    % set color for negative values
    nendColor = [1 1 1]; % should not be completely white, otherwise hits with low S/N are almost not visible
    nstartColor = [0.9 0.9 0.9];
    
    if  threshold
        % determine location of threshold
        if threshold < minval
            mintpos = round(numOfGradations*threshold/(minval+maxval),0);
            negcont = neggrad - mintpos; 
        else
            mintpos = minval;
            negcont = 0;
        end
        if threshold < maxval
            maxtpos = round(numOfGradations*threshold/(minval+maxval),0);
            poscont = posgrad -maxtpos;
        else
            maxtpos = minval;
            poscont = 0;
        end

        negcont = get_gradient(nstartColor,nstartColor,negcont);
        neggradient = get_gradient(nstartColor,nendColor,mintpos);
        posgradient =  get_gradient(startColor,endColor,maxtpos);
        poscont = get_gradient(endColor,endColor,poscont);
        colormap([negcont;neggradient;1 1 1; posgradient;poscont]);  % Add white color at 0 - for when SN=0;
 
    else
        neggradient = get_gradient(nstartColor,nendColor,neggrad);
        posgradient = get_gradient(startColor,endColor,posgrad);
        colormap([neggradient;1 1 1; posgradient]);  % Add white color at 0 - for when SN=0;
    end
else
    if  threshold
        
        % for area below 0
        nstartColor = [1 1 1];
        negcont = minval;
        
        % determine location of threshold
        if threshold < maxval
            maxtpos = round(numOfGradations*threshold/(minval+maxval),0);
            poscont = posgrad -maxtpos;
        else
            maxtpos = minval;
            poscont = 0;
        end
        
        negcont = get_gradient(nstartColor,nstartColor,negcont);
        posgradient =  get_gradient(startColor,endColor,maxtpos);
        poscont = get_gradient(endColor,endColor,poscont);
        colormap([negcont;1 1 1; posgradient;poscont]);  % Add white color at 0 - for when SN=0;

    else
        gradient = get_gradient(startColor,endColor,numOfGradations);
        colormap([1 1 1; gradient]);  % Add white color at the start - for when SN=0;
        
    end
end


end


function grad = get_gradient(startColor,endColor,steps)

grad = [linspace(startColor(1),endColor(1),steps);...
    linspace(startColor(2),endColor(2),steps);...
    linspace(startColor(3),endColor(3),steps)]';

end

function [yCalc,R2] = fitlnmodel(x,y)

X = [ones(length(x),1) x]; % add 1 for y-intercept
p = X\y;

yCalc = X*p;
R2 = 1 - sum((y - yCalc).^2)/sum((y - mean(y)).^2);

end

function plot_overview_interactions(binmapfull,dims,Names,results_folder,nametag)

if find(strcmp(Names, 'ATP'))
    load('input_other\mNames_sorted.mat');
    pNames = mNames_sorted;
    [~,~,pos2] = intersect(pNames,Names,'stable');
    sortmap = sum(binmapfull,dims);
    sortmap = sortmap(pos2);
else
    [sortmap,sortidx] = sort(sum(binmapfull,dims));
    pNames = Names(sortidx);
end
f = figure('visible','off');
set(f, 'Color', [1,1,1], 'Position', [100 100 length(Names)*100 500]); % [x y width height]
hold all;
set(gca,'FontSize', 14,'LineWidth',2);
h = bar(sortmap,'LineWidth',2);

set(gca,'xTick',[1:length(sortmap)]);
set(gca,'xticklabel',pNames);
set(gca,'tickdir','out');
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 14)
xtickangle(90);


ylabel('# of detected interactions');
title(char(strcat('Detected interactions per ',nametag)));
hold off;
fname = fullfile(results_folder,char(strcat('Interactions_',nametag,'_bars.png')));
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(f,fname);

% also as pie chart
pNames(sortmap==0) = [];
sortmap(sortmap==0) = [];
for i = 1:length(pNames)
    if sortmap(i) > 1
        pNames{i} = char(strcat(pNames{i},' (',num2str(sortmap(i)),')'));
    end
end
f = figure('visible','off');
set(f, 'Color', [1,1,1], 'Position', [0 100 600 600]); % [x y width height]
hold all;
set(gca,'FontSize', 10);

p = pie(sortmap,pNames);
ax = gca;
ax.Visible = 'off';
dists = [0.91, 0.96];
for iHandle = 2:2:2*numel(pNames)
    if sortmap(iHandle/2) == 1
        p(iHandle).Position = dists(mod(iHandle/2,2)+1)*p(iHandle).Position;
    else
        p(iHandle).Position = 0.92*p(iHandle).Position;
    end
end
title(['Detected interactions per ' nametag]);
hold off;
fname = fullfile(results_folder,char(strcat('Interactions_',nametag,'_pie.png')));
saveas(f,fname);

end

function plot_percentual_interactions(binmapfull,mNames, pNames, results_folder,knownInteractions)

% check for every how much is known
reg = struct;
reg.met = zeros(length(mNames),2);
reg.prot = zeros(length(pNames),2);

cat = struct;
cat.met = zeros(length(mNames),2);
cat.prot = zeros(length(pNames),2);


for i = 1:length(mNames) 
    for j = 1:length(pNames)
 
        pidx = find(strcmpi(knownInteractions.prot,pNames(j)));
        midx = find(strcmpi(knownInteractions.met(pidx,:),mNames(i)));
        
        if length(midx) > 0
            
            if strfind(knownInteractions.type{pidx,midx},'S')
                cat.met(i,1) = cat.met(i,1) + 1;
                cat.prot(j,1) = cat.prot(j,1) + 1;
                
                if binmapfull(j,i) > 0
                    cat.met(i,2) = cat.met(i,2) + 1;
                    cat.prot(j,2) = cat.prot(j,2) + 1;
                end
            end
            
            if strfind(knownInteractions.type{pidx,midx},'R')
                reg.met(i,1) = reg.met(i,1) + 1;
                reg.prot(j,1) = reg.prot(j,1) + 1;
                
                if binmapfull(j,i) > 0
                    reg.met(i,2) = reg.met(i,2) + 1;
                    reg.prot(j,2) = reg.prot(j,2) + 1;
                end
            end
        end
        
    end
end

% compute percentages
reg.met(:,3) = reg.met(:,1) - reg.met(:,2);
cat.met(:,3) = cat.met(:,1) - cat.met(:,2);
reg.prot(:,3) = reg.prot(:,1) - reg.prot(:,2);
cat.prot(:,3) = cat.prot(:,1) - cat.prot(:,2);

plot_bars(reg.met, mNames, results_folder, 'regulatory', 'metabolite');
plot_bars(reg.prot, pNames,results_folder, 'regulatory', 'protein');
plot_bars(cat.met, mNames, results_folder, 'catalytic', 'metabolite');
plot_bars(cat.prot, pNames, results_folder, 'catalytic', 'protein');

end


function plot_bars(matrix, Names, results_folder, nametag1, nametag2)

if find(strcmp(Names, 'ATP'))
    load('input_other\mNames_sorted.mat');
    pNames = mNames_sorted;
    [~,~,pos2] = intersect(pNames,Names,'stable');
    matrix = matrix(pos2,:);
    Names = pNames;
end

%keeps = find(matrix(:,1) > 0);
%Names = Names(keeps);
%matrix = matrix(keeps,:);
cmaps = [0.22 0.62 0; 1 1 1];

percs(:,1) = matrix(:,2) ./ matrix(:,1); 
percs(:,2) = matrix(:,3) ./ matrix(:,1); 

f = figure('visible','off');
set(f, 'Color', [1,1,1], 'Position', [100 100 length(Names)*100 500]); % [x y width height]
hold all;

set(gca,'FontSize', 14,'LineWidth',2);

h = bar(matrix(:,2:3), 'stacked', 'FaceColor', 'flat','LineWidth',2);
for k = 1:size(h,2)
    h(k).CData = cmaps(k,:);
end
%legend(h,{'detected', 'not detected'}, 'location', 'northwest');

set(gca,'xTick',[1:length(Names)]);
set(gca,'xticklabel',Names);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 14)
xtickangle(90);
set(gca,'tickdir','out');

ylim([0 8]);

ylabel('# detected interactions');
title(char(strcat({'Known, '}, nametag1, {' interactions per '},  nametag2)));
hold off;
fname = fullfile(results_folder,char(strcat('Bars_abs_known_',nametag1,'_',nametag2,'.png')));
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
saveas(f,fname);


f = figure('visible','off');
hold all;

set(gca,'FontSize', 12);
h = bar(percs, 'stacked', 'FaceColor', 'flat');
for k = 1:size(h,2)
    h(k).CData = cmaps(k,:);
end

set(gca,'xTick',[1:length(Names)]);
set(gca,'xticklabel',Names);
xAX = get(gca,'XAxis');
set(xAX,'FontSize', 10)
xtickangle(90);

ylabel('% detected interactions');
title(char(strcat({'Known, '}, nametag1, {' interactions per '},  nametag2)));
hold off;
fname = fullfile(results_folder,char(strcat('Bars_perc_known_',nametag1,'_',nametag2,'.png')));
saveas(f,fname);


end

function DisplayHitMap(absolute_colormap_threshold,fSize,mapCellSize,zscores,mapFull,pNames,mNames,knownInteractions,stats,fnamebase,mixes)

n_compounds = numel(mNames);
n_prot = numel(pNames);

% change order if mix info is there
if nargin > 10
   [~,sortidx] = sort(mixes);
   mNames = mNames(sortidx);
   mapFull = mapFull(:,sortidx);
else
    load('input_other\mNames_sorted.mat');
    [~,~,pos2] = intersect(mNames_sorted,mNames,'stable');
    mapFull = mapFull(:,pos2);
    mNames = mNames_sorted;
end


f = figure('visible','off');
rows = 2; % need to expand figue size - to fit metabolite names
set(f, 'Color', [1,1,1], 'Position', [0 100 n_compounds*mapCellSize n_prot*mapCellSize*rows]); % [x y width height]

%%% compute colormap
compute_colormap(zscores,mapFull,absolute_colormap_threshold);

plot_handle = subplot(4,5,[7:10,12:15,17:20]); % need to expand figue size - to fit metabolite names
POS = get(plot_handle, 'pos');
imagesc(mapFull);

vertLabels = pNames;
set(gca,'YTick',[1:length(pNames)]);
set(gca,...
    'XAxisLocation','top',...
    'TickDir', 'in',...
    'Ticklength', [0 0],...
    'YTickLabel', vertLabels,...
    'XTickLabel', [],... % could use compoundNames here
    'XTickLabelMode', 'manual',...
    'FontSize', fSize...
    );

colorbar('eastoutside', 'FontSize', fSize-1);


    % plot * on known inetractions
     for i=1:length(pNames)
        pidx = find(strcmpi(knownInteractions.prot,pNames(i)));
        for j = 1:length(mNames)
            midx = find(strcmpi(knownInteractions.met(pidx,:),mNames(j)));
            if length(midx) > 0
                % plot *
                txt ='';
                tsize = 17;
                if strfind(knownInteractions.type{pidx,midx},'S')
                    txt = strcat(txt,'C');
                    tsize = tsize-5;
                end
                if strfind(knownInteractions.type{pidx,midx},'R')
                    txt = strcat(txt,'R');
                    tsize = tsize-5;
                end
                
                if mapFull(i,j) > 0
                    text(j-0.3,i,txt,'FontSize',tsize,'Color','black');
                else
                    text(j-0.3,i,txt,'FontSize',tsize,'Color','red');
                end
                
            end
        end
     end



% plot lines if mixes are plotted
if nargin > 10
    lengths = [];
    for i = 1:3
       lengths(i) = length(find(mixes==i)); 
       v1(1:length(pNames)+2) = sum(lengths(1:i))+0.5;
       line(v1,0:length(pNames)+1,'Color','red','LineWidth',1);
    end
end


% text(X,Y,string);
textY = POS(2)+POS(4)-0.3;
% text(1:n_compounds,repmat( mapCellSize*0.01,1,n_compounds), mNames,'HorizontalAlignment','left','FontSize',fSize,'rotation',90);
text(1:n_compounds,repmat( textY,1,n_compounds), mNames,'HorizontalAlignment','left','FontSize',fSize,'rotation',90);

% if Iavailable -> plot psizes

    plot_handle = subplot(4,5,[6,11,16]);
    grad = get_gradient([0.33 0.4 0.96],[0.85 0.86 0.97],8);
    %colormap(plot_handle,[grad]);
    d = flipud(stats.protSize); % flip to have the same orientation as in matrix
    %h = barh(d,'stacked'); 
    h = barh(d,'stacked','FaceColor','flat');
    
    %color
    for k = 1:size(d,2)
        h(k).CData = grad(k,:);
    end
    
   % plot line at 50 kDa
    v1(1:length(pNames)+2) = 50;
    line(v1,0:length(pNames)+1,'Color','red','LineWidth',1);
    
    set(gca,'ylim',[0.5 length(pNames)+0.5]);
    set(gca,'yTick',[1:length(pNames)]);
    set(gca,'yticklabel',fliplr(pNames));
    xmax = max(sum(d'));
    set(gca,'xlim',[0 xmax+10]);
    xlabel('kDa');
    title('Protein Size');
    %breakxaxis([225 475],0.005);   


%%% Save hit maps as pdfs
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
fName = char(strcat(fnamebase,'.pdf'));
saveas(f,fName);
%print(f,fnamebase,'-dpdf','-r0');

%%% also save as png
fName = char(strcat(fnamebase,'.png'));
saveas(f,fName);

end