function Generate_Interaction_Map_MakeMetaboliteOverview(mapFull,pNames,mNames,knownInteractions,results_folder)

%change results folder to subfolder
results_folder = char(strcat(results_folder,'/overview_detected_interactions/'));
if ~exist(results_folder, 'dir')
    mkdir(results_folder); % folder for output
    disp('Created folder for overview output');
end


compounds(1).names = {'ATP','GTP','CTP','TTP','AMP', 'IMP', 'cAMP'};
compounds(2).names = {'ALA','ARG','ASN','ASP','CYSox','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'};
compounds(3).names = {'3M2O','CYSTAT','HSER','HYXA','MALON','ORN','PRPP','SKM','SPER','URA','ITA','PhePYR'};
compounds(4).names = {'FAD','NAD','NADP'};
compounds(5).names = {'3PG','6PG','AKG','CIT','DHAP','FBP','FUM','G6P','MAL','PEP','PYR','R5P','SUC'};


groupname = {'nucleotides','aminoacids','others','cofactors','CCM'};

txt = 0; %1:plot text, 2: plot dots

% compute matrix of known interactions
knownmap = zeros(length(pNames),length(mNames));
for j = 1:length(pNames)
    pidx = find(strcmpi(knownInteractions.prot, pNames(j)));
    for i = 1:length(mNames)
        midx = find(strcmpi(knownInteractions.met(pidx,:),mNames(i)));
         if ~isempty(midx)
             if strfind(knownInteractions.type{pidx,midx},'S')
                 knownmap(j,i) = 2;
             end
             if strfind(knownInteractions.type{pidx,midx},'R')
                 knownmap(j,i) = 1;
             end
         end
    end
end


fSize = 12;

%%% Display table concerning all cpds
for i = 1:length(compounds)
    % extract submap
    [~,pos1,pos2] = intersect(mNames,compounds(i).names);
    new_mNames = compounds(i).names(pos2);
    new_map = mapFull(:,pos1);
    new_knownmap = knownmap(:,pos1);
    
    kpidx = find(or(sum(new_map,2)>0,sum(new_knownmap,2)>0));
    new_map = new_map(kpidx,:);
    new_pNames = pNames(kpidx);
    kpidx = find(or(sum(new_map,1)>0,sum(new_knownmap,1)>0));
    new_map = new_map(:,kpidx);
    new_mNames = new_mNames(kpidx);
    
    filename = char(strcat(results_folder,'/interactions_with_',groupname(i),'.pdf'));
    %disp_list(new_mNames, new_pNames, new_map, fSize, knownInteractions, txt, filename);    
    disp_list(new_mNames, new_pNames, new_map, fSize, new_knownmap, txt, filename);    
end




end

function disp_list(new_mNames, new_pNames, new_map, fSize,knownmap, txt, fName)


n_compounds = numel(new_mNames);
n_prot = numel(new_pNames);

% compute colormatrix
if n_compounds < 8
   colormatrix = colormap(lines(n_compounds)); 
else
   colormatrix = colormap(parula(round(n_compounds*1.3,0)));
end

% distance for plotting
dists(1) = 0.5;
dists(2) = dists(1)+0.5;

% make figure
f = figure('visible','off');
set(f, 'Color', [1,1,1], 'Position', [100 100 (n_prot+1)*120 (n_compounds)*100]); % [x y width height]
hold all;
scatter([0 size(new_map,1)],[0 size(new_map,2)]);

for j = 1:n_prot
    
    % write protein name
    text(j-0.3,size(new_map,2),new_pNames(j),'FontSize',fSize,'Color','black')
    dctr = dists(mod(j,2)+1);
    hctr = -0.5;
    
    for i = 1:n_compounds
        
        % find if interaction is known
        currstr = new_mNames(i);
        knownflag = 0;
        knownflagr = 0;
        
        if knownmap(j,i) == 1 % regulator
            currstr = strcat(currstr,'*R');
            knownflag = 1;
            knownflagr = 1;
        elseif knownmap(j,i) == 2 % substrate
            currstr = strcat(currstr,'*S');
            knownflag = 1;
        end
        
        % check if interaction was detected and print
        if new_map(j,i) > 0
            if txt
                text(j-0.3,size(new_map,2)-dctr, currstr,'FontSize',fSize,'Color',colormatrix(i,:));
                dctr = dctr +dists(1);
            else
                if knownflag
                    plot(j+hctr,size(new_map,2)-dctr, 'o', 'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',colormatrix(i,:),'MarkerSize',12, 'LineWidth', 3);
                else
                    plot(j+hctr,size(new_map,2)-dctr, 'o', 'MarkerFaceColor',colormatrix(i,:),'MarkerEdgeColor',colormatrix(i,:),'MarkerSize',12, 'LineWidth', 3);
                end
                hctr = hctr+0.28;
            end
        else
            if knownflagr
                plot(j+hctr,size(new_map,2)-dctr, 'o', 'MarkerFaceColor',[1 1 1],'MarkerEdgeColor',colormatrix(i,:),'MarkerSize',12, 'LineWidth', 3);
                hctr = hctr+0.28;
            end
        end
    end
end

% plot legend
if ~txt
    dctr = dists(2)*2 + (dists(1));
    for i = 1:n_compounds
       text(0.7,dctr, new_mNames(i),'FontSize',fSize,'Color',colormatrix(i,:));
       plot(2,dctr, '.', 'Color',colormatrix(i,:),'MarkerSize',40);
       dctr = dctr -(dists(1)/2);
    end
end

ax = gca;
ax.Visible = 'off';

%%% Save hit maps as pdfs
set(f,'Units','Inches');
pos = get(f,'Position');
set(f,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)]);
print(f,fName,'-dpdf','-r0');

end
