function Stats_AnalyzeProteinStability

% function to calculate the stablity of the proteins during NMR measurements
%
% the script calculates the integral of specified regions in the NMR spectrum
% ('aromatic' and 'aliphatic' regions),
% finally summing up all the values. If this integral is constant over 
% multiple experiments with the same protein, the protein is
% considered to be stable
%
%% !!!
% this script requires the 'rbnmr' script 
% https://ch.mathworks.com/matlabcentral/fileexchange/40332-rbnmr)
% 
% input
% the script expects an xlsx wich indicates
%       on sheet 1 ('ppm regions'): the ppm regions to investigate
%       on sheet 2 ('datasets'):    where the spectra of 'pure proteins' can be found ('real' spectra, not exported txt files)
% see 'input_other/PPM_ProteinStability.xlsx'
%
% output
% results figures are saved in 'resfolder'

%% set params
% input
proteinfile = 'input_other/PPM_ProteinStability.xlsx';

% output
resfolder = 'results_proteinstability';

if ~exist(resfolder, 'dir')
    mkdir(resfolder); % folder for output
    disp('Created folder for output');
end
%%

% read region information
region = struct;
[region.left,region.type,~] = xlsread(proteinfile,'ppm regions','A1:E500');
region.right = region.left(:,2);
region.left(:,2) = []; 

% read dataset information
ds = dataset_information(proteinfile,'datasets');

for i = 1:length(ds)
    ds(i).spectrum = ReadSpecFile(ds(i).path,10,-0.5);
    ds(i).integral = CalculateIntegral(ds(i).spectrum,region);
end

% transform to datamatrix
[matrix,mixnames,pnames] = ToMatrix(ds);

% plot bar graph
f = figure('visible','off');
c = categorical(pnames);
bar(c,matrix');
legend(mixnames);
fName = fullfile(resfolder,'ProteinStability.png');
saveas(f,fName);

% compute stds
smatrix = std(matrix) ./ mean(matrix);
f = figure('visible','off');
bar(c,smatrix');
ylim([0 1]);
title('Std(all mixes) / mean(all mixes)');
fName = fullfile(resfolder,'ProteinStability_Std.png');
saveas(f,fName);

end


function [matrix,mixnames,pnames] = ToMatrix(ds)

matrix = reshape([ds.integral],5,[]);
pnames = unique([ds.protein],'stable');
mixnames = unique([ds.mix],'stable');

end


function ds = dataset_information(proteinfile,sheet)

ds = struct;
mixes = {'nomix','ccm1','ccm2','ccm3','ccm4'};
mixidx = [2,12,22,32,42];

[dsnums,datasets,~] = xlsread(proteinfile,sheet,'A1:E500');

ctr = 1;
for i = 1:length(dsnums)
    for j = 1:length(mixes)
        ds(ctr).protein = datasets(i,1);
        ds(ctr).mix =mixes(j);
        ds(ctr).dsnum = dsnums(i)+mixidx(j); 
        ds(ctr).path = fullfile(datasets(i,3),datasets(i,4),num2str(ds(ctr).dsnum),'\pdata\1');
        ctr = ctr+1;
    end
end


end

function res=ReadSpecFile(path,xleft,xright)

% function that reads NMR files
% input
% path: path to Spec file
% xleft: leftmost x-value to include
% xright: rightost x-value to include
% output
% res: structre with fields: info, x, data, data_rescaled, xleft, xright

%1: read file
tmp = rbnmr(path);

res.info = tmp.Info.FilePath;

res.xleft = xleft;
res.xright = xright;

% find pos of xleft and xright
tmpx = tmp.XAxis - xleft;
[~,xindleft] = min(abs(tmpx));
tmpx = tmp.XAxis - xright;
[~,xindright] = min(abs(tmpx));

res.data=tmp.Data(xindleft:xindright);
res.xaxis=tmp.XAxis(xindleft:xindright);
res.data_rescaled=res.data./max(res.data);

end


function integral = CalculateIntegral(spec,region)

ydata = spec.data;
xdata = spec.xaxis;
integral = 0;

for i = 1:length(region.left)

    % find left & rightx
    [~,lx] = min(abs(xdata-region.left(i)));
    [~,rx] = min(abs(xdata-region.right(i)));
    
    % calculate integral
    cint = trapz(xdata(lx:rx),ydata(lx:rx));
    if strcmp(region.type{i},'Aromatic')
        integral = integral - cint;
    else
        integral = integral + cint;
    end
    
end

end