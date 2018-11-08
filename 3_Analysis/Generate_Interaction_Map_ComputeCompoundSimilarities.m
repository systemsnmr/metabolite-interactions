function Generate_Interaction_Map_ComputeCompoundSimilarities(mapFull,cpdsim,mode,pNames,mNames,knownInteractions,results_folder)

% change metabolite name order
load('input_other\mNames_sorted.mat');
[~,~,pos2] = intersect(mNames_sorted,mNames,'stable');
mapFull = mapFull(:,pos2);
mNames = mNames_sorted;

%change results folder to subfolder
results_folder = char(strcat(results_folder,'/compound_similarities_',mode,'/'));
if ~exist(results_folder, 'dir')
    mkdir(results_folder); % folder for output
    disp('Created folder for similarity output');
end

fnames = [cpdsim.info.name];

sim  = struct;
sim.hits.cumsims = [];
sim.hits.maxsims = [];
sim.hits.simmat = zeros(length(pNames),length(mNames)); % matrix of max similarities

sim.all.maxsims = [];
sim.all.simmat = zeros(length(pNames),length(mNames)); % matrix of all max similarities

sim.not.maxsims = []; % for similarity of known but not detected
sim.newhits.maxsims = []; % for similarity of new, not known

sim.prot.regknown = 0; % for known regulatory inetractions per protein
sim.prot.new = 0; % for number of new interactions per protein
sim.prot.newlowsim = 0; % for number of new interactions per protein with low similarity
simcutoff = 0.5; % cutoff for newlosim

% metabolites in mixes
midx = find_idx(cpdsim.metmixes(:,1), mNames); 

% iterate over all proteins
for i = 1:length(pNames)
    
    sim.prot(i).regknown = 0;
    sim.prot(i).new = 0; 
    sim.prot(i).newlowsim = 0;
    
    % find substrates indices
    pidx = find(strcmpi(fnames,pNames{i}));
    substrates = cpdsim.info(pidx).catalytic;
    sidx = find_idx(cpdsim.allmets(:,1), substrates);
    
    allsims = [];

    % check for known but not detected interactions
    ppidx = find(strcmpi(knownInteractions.prot,pNames(i)));
    
    for j = 1:size(mapFull,2)
        
        % find similarities with susbtrates for every metabolite
%        for k = 1:length(midx)
         allsims(j,:) = cpdsim.(mode).simmat(sidx,midx(j));  
%        end
         mmidx = find(strcmpi(knownInteractions.met(ppidx,:),mNames(j)));
         if mapFull(i,j) == 0
            if length(mmidx) > 0 % known, but not detected
                sim.not.maxsims = [sim.not.maxsims max(allsims(j,:),[],2)'];
            end
         else
            if length(mmidx) == 0 % not known, but detected
                sim.newhits.maxsims = [sim.newhits.maxsims max(allsims(j,:),[],2)'];
            end
         end
         
         % count!
         if length(mmidx) > 0 % known
             if strfind(knownInteractions.type{ppidx,mmidx},'R')
                 sim.prot(i).regknown = sim.prot(i).regknown + 1; % all knowns
             end
         else % unknown
             if mapFull(i,j) > 0
                 sim.prot(i).new = sim.prot(i).new + 1; % all new
                 if max(allsims(j,:)) < simcutoff
                     sim.prot(i).newlowsim = sim.prot(i).newlowsim +1; %new & low sim
                 end
             end
         end
        
    end
    
    % save info on all possible
    sim.all.simmat(i,:) = max(allsims,[],2);
    sim.all.maxsims = [sim.all.maxsims max(allsims,[],2)'];
        
    % find all metabolites with score > 0
    idx = find(mapFull(i,:) > 0);
    
    if ~isempty(idx)
        relsims =  allsims(idx,:);
        
        % save info on hits
        sim.hits.cumsims = [sim.hits.cumsims reshape(relsims,1,[])];
        sim.hits.maxsims = [sim.hits.maxsims max(relsims,[],2)'];
        sim.hits.simmat(i,idx) = sim.all.simmat(i,idx);
        
        filename = char(strcat(results_folder,'/',mode,'_similarities_for_',pNames(i),'.png'));
        plot_similarities(pNames(i), mNames(idx), relsims, filename);
    end
    
end


% plot histograms
cmat = [0,0,1;1,0,0;0,1,0];
names = {'all possible','detected'};
m1 = reshape(cpdsim.(mode).simmat,1,[]);
m2 = sim.hits.cumsims;

titlen = 'Similarity of all possible interactions vs. detected interactions';
filename = char(strcat(results_folder,'/',mode,'_similarities_overview.png'));
plot_histograms(cmat,names,0,m1,m2,titlen,filename);


titlen = 'Normalized similarity of all possible interactions vs. detected interactions';
filename = char(strcat(results_folder,'/',mode,'_similarities_overview_norm.png'));
plot_histograms(cmat,names,1,m1,m2,titlen,filename);


% max detected vs not detected
cmat = [0,0,1;0,1,0];
names = {'max(all possible)','max(detected)'};
m1 = sim.all.maxsims;
m2 = sim.hits.maxsims;

titlen = 'Max similarity of all possible interactions vs. detected interactions';
filename = char(strcat(results_folder,'/',mode,'_similarities_max.png'));
plot_histograms(cmat,names,0,m1,m2,titlen,filename);

titlen = 'Normalized max similarity of all possible interactions vs. detected interactions';
filename = char(strcat(results_folder,'/',mode,'_similarities_max_norm.png'));
plot_histograms(cmat,names,1,m1,m2,titlen,filename);

% plot map with similarity scores
filename = char(strcat(results_folder,'/',mode,'_overview'));
plot_map(pNames,mNames,sim.hits.simmat,knownInteractions,filename,mapFull);
% save also as xlsx
filename = char(strcat(results_folder,'/',mode,'_table.xlsx'));
xlswrite(filename,sim.hits.simmat);

filename = char(strcat(results_folder,'/',mode,'_overview_all'));
plot_map(pNames,mNames,sim.all.simmat,knownInteractions,filename);

% plot with not detected
cmat = [0,0,1;0,1,0;1,0,0];
names = {'max(all possible)','max(detected)','max(known,not detected)'};
m3 = sim.not.maxsims;
titlen = 'Normalized max similarity';
filename = char(strcat(results_folder,'/',mode,'_similarities_max_norm_addinfo.png'));
plot_histograms(cmat,names,1,m1,m2,titlen,filename,m3);


% save distribution of similarities of new hits as xlsx
filename = char(strcat(results_folder,'/',mode,'_similarities_max_newhits.xlsx'));
xlswrite(filename,sort(sim.newhits.maxsims)');

m3 = sim.newhits.maxsims;
titlen = 'Normalized max similarity of all possible interactions vs. detected interactions and news';
filename = char(strcat(results_folder,'/',mode,'_similarities_max_news_norm.png'));
plot_histograms(cmat,names,1,m1,m2,titlen,filename,m3);

% save sim in struct
filename = char(strcat(results_folder,'/',mode,'_sim.mat'));
save(filename,'sim');


%% plot of count new hits and known hits

sortmap = [[sim.prot.regknown]; [sim.prot.new]; [sim.prot.newlowsim]]';

cmat = [0,0,1;0,1,0;1,0,0];
names = {'known reg interactions','new interactions','new interactions, sim < 0.5'};

f = figure('visible','off');
hold all;
set(gca,'FontSize', 12);
h = barh(sortmap);

set(gca,'yTick',[1:length(sortmap)]);
set(gca,'yticklabel',pNames);
yAX = get(gca,'YAxis');
set(yAX,'FontSize', 10)

legend([h],[names]);

hold off;
fname = fullfile(results_folder,char(strcat('Bars_interaction.png')));
saveas(f,fname);

end


function plot_histograms(cmat,names,norm,m1,m2,titlen,filename,m3)

f = figure('visible','off');
hold all;
set(gca,'FontSize', 14,'LineWidth',2);

if not(norm)
    a(1) = histogram(m1,'FaceColor',cmat(1,:),'LineWidth',2);
    a(2) = histogram(m2,'FaceColor',cmat(2,:),'LineWidth',2);        
    a(1).BinWidth = 0.025;
    a(2).BinWidth = 0.025;
    
    if nargin > 7
        a(3) = histogram(m3,'FaceColor',cmat(3,:),'LineWidth',2);        
        a(3).BinWidth = 0.025; 
    end
    
    ylim([0 150]);
    znums = length(find(m1 ==0));
    text(0.01,150,char(strcat('Max:',num2str(znums))));
    ylabel('counts');
else
    a(1) = histogram(m1,'FaceColor',cmat(1,:),'normalization','probability','LineWidth',2);
    a(2) = histogram(m2,'FaceColor',cmat(2,:),'normalization','probability','LineWidth',2);
    a(1).BinWidth = 0.025;
    a(2).BinWidth = 0.025;
    
    if nargin > 7
         a(3) = histogram(m3,'FaceColor',cmat(3,:),'normalization','probability','LineWidth',2);     
         a(3).BinWidth = 0.025; 
    end
    
    ylabel('%');
end

legend([a],[names]);
set(gca,'tickdir','out');
xlabel('similarity');
title(titlen);
hold off;
saveas(f,filename)
end



function plot_similarities(pName, mNames, allsims, filename)

% assemble mat for scattering
scatmat = [];
for i = 1:size(allsims,1)
    for j = 1:size(allsims,2)
        scatmat(end+1,1) = i;
        scatmat(end,2) = allsims(i,j);
    end
end

% boxplot of similarities
f = figure('visible','off');
hold all;
boxplot(allsims','Labels',mNames);
scatter(scatmat(:,1),scatmat(:,2),'filled','k');
title(char(strcat('Similarities of regulators with substrates of:',pName)));
ylim([0 1]);
hold off;
saveas(f,filename);
end

function idx = find_idx(allmets, cpds)
idx = [];
for i = 1:length(cpds)
    idx(i) = find(strcmp(allmets,cpds(i)));
end
end


function plot_map(pNames,mNames,mapFull,knownInteractions,fnamebase,fsimap)

n_compounds = numel(mNames);
n_prot = numel(pNames);
fSize = 11;

%%% Display hitMap
f = figure('visible','off');
rows = 2; % need to expand figue size - to fit metabolite names
set(f, 'Color', [1,1,1], 'Position', [0 100 n_compounds*30 n_prot*30*rows]); % [x y width height]

%%% compute colormap
gradient = get_gradient([0.9 0.9 0.95],[0.35 0.35 0.95],1000);
colormap([1 1 1; gradient]);  % Add white color at the start - for when =0;

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

simplot = 1;
% plot * on known inetractions or by similarity
if simplot == 1
    tsize = 17;
    if nargin > 5
        for i=1:length(pNames)
            for j = 1:length(mNames)
                if and(mapFull(i,j) <= 0.5,fsimap(i,j) > 0)
                    rectangle('Position',[j-0.5 i-0.5 1 1],'EdgeColor','black','LineWidth',2);
                    %text(j-0.3,i,'x','FontSize',tsize,'Color','black');
                end
            end
        end
    else
       for i=1:length(pNames)
            for j = 1:length(mNames)
                if mapFull(i,j) <= 0.5
                    rectangle('Position',[j-0.5 i-0.5 1 1],'EdgeColor','black','LineWidth',2);
                    %text(j-0.3,i,'x','FontSize',tsize,'Color','black');
                end
            end
        end
    end 
        
else
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
end

% text(X,Y,string);
textY = POS(2)+POS(4)-0.3;
% text(1:n_compounds,repmat( mapCellSize*0.01,1,n_compounds), mNames,'HorizontalAlignment','left','FontSize',fSize,'rotation',90);
text(1:n_compounds,repmat( textY,1,n_compounds), mNames,'HorizontalAlignment','left','FontSize',fSize,'rotation',90);

%%% Save hit maps as pdf and png
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f,fnamebase,'-dpdf','-r0');

fName = char(strcat(fnamebase,'.png'));
saveas(f,fName);


end

function grad = get_gradient(startColor,endColor,steps)

grad = [linspace(startColor(1),endColor(1),steps);...
    linspace(startColor(2),endColor(2),steps);...
    linspace(startColor(3),endColor(3),steps)]';

end