function Generate_Interaction_Map_MakeOutliersFigures(fSize,mapCellSize,addinfo,pNames,mNames,stats,results_folder)

% make folder for results
results_folder = char(strcat(results_folder,'\Neg_S1S2_Outliers\'));
if ~exist(results_folder, 'dir')
    mkdir(results_folder);
end

fName = fullfile(results_folder,'Map_NegFlag_');
plot_OutliersMap(addinfo.hNeg,pNames,mNames,addinfo.peakNames,fSize,mapCellSize,stats, fName);

fName = fullfile(results_folder,'Map_S1S2Flag_');
plot_OutliersMap(addinfo.hS1S2,pNames,mNames,addinfo.peakNames,fSize,mapCellSize,stats, fName);


end

function plot_OutliersMap(currstruct,pNames,mNames,peakNames,fSize,mapCellSize,stats,namebase)
    
n_prot = numel(pNames);
mix = unique({currstruct.mix});

 for i = 1:length(mix)
       currh = currstruct(strcmp({currstruct.mix},mix(i))); 
       currpeaks = peakNames(i).Names;%unique({currh.peakid});
       n_compounds = numel(currpeaks);
       currmap = zeros(length(pNames),length(currpeaks)); 
       fName = char(strcat(namebase,mix(i),'.png'));
       
       % first, compute interaction map
       for pIdx = 1:length(pNames)
           for mIdx = 1:length(currpeaks)
               % find all matches S1s2
               allpids = strcmp({currstruct.protein},pNames(pIdx));
               allmids = strcmp({currstruct.peakid},currpeaks(mIdx));
               allids = find((allpids+allmids)==2);
               if ~isempty(allids)
                   currmap(pIdx,mIdx) = 1;
               end
           end
       end

       %%% Display hitMap
      Display_Map(n_compounds,fSize, mapCellSize,n_prot,currmap,pNames,currpeaks,stats,fName);
 end
 
 currpeaks = mNames;
 n_compounds = numel(currpeaks);
 currmap = zeros(length(pNames),length(currpeaks));
 fName = char(strcat(namebase,'.png'));
 
 % first, compute interaction map
 for pIdx = 1:length(pNames)
     for mIdx = 1:length(currpeaks)
         % find all matches S1s2
         allpids = strcmp({currstruct.protein},pNames(pIdx));
         allmids = strcmp({currstruct.Assign},currpeaks(mIdx));
         allids = find((allpids+allmids)==2);
         if ~isempty(allids)
             currmap(pIdx,mIdx) = 1;
         end
     end
 end
 
 %%% Display hitMap
 Display_Map(n_compounds,fSize, mapCellSize,n_prot,currmap,pNames,currpeaks,stats,fName);
 
 
end

function compute_colormap
startColor = [0.84 1 0.76]; % should not be completely white, otherwise hits with low S/N are almost not visible
endColor = [0.22 0.62 0];
steps = 1000;

gradient = [linspace(startColor(1),endColor(1),steps);...
    linspace(startColor(2),endColor(2),steps);...
    linspace(startColor(3),endColor(3),steps)]';


colormap([1 1 1; gradient]);  % Add white color at the start - for when SN=0;
       
end

function grad = get_gradient(startColor,endColor,steps)

grad = [linspace(startColor(1),endColor(1),steps);...
    linspace(startColor(2),endColor(2),steps);...
    linspace(startColor(3),endColor(3),steps)]';

end

function Display_Map(n_compounds,fSize, mapCellSize,n_prot,currmap,pNames,currpeaks,stats,fName)

f = figure('visible','off');
rows = 2; % need to expand figue size - to fit metabolite names
set(f, 'Color', [1,1,1], 'Position', [0 100 n_compounds*mapCellSize n_prot*mapCellSize*rows]); % [x y width height]

%%% compute colormap
compute_colormap;

plot_handle = subplot(4,5,[7:10,12:15,17:20]); % need to expand figue size - to fit metabolite names
POS = get(plot_handle, 'pos');
imagesc(currmap);

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

textY = POS(2)+POS(4)-0.3;
text(1:n_compounds,repmat( textY,1,n_compounds), currpeaks,'HorizontalAlignment','left','FontSize',fSize,'rotation',90);

% -> plot psizes
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

set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);

%%% save as png
saveas(f,fName);
end