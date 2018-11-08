function out =  Stats_ParameterScanROC

% script to scan the influence of several parameters on the
% false-positive/true positive rate of recovered interactions

%%
% define parameter range to scan
fsamode =1:1:2; % 1: fractional signal intensity (FSI)
                % 2: relaxation factor (RF)
scoringopt = 1; % 1: use s/n cutoff and lower cutoff
                % 2: use rank based scoring (equal weight on S/N and FSA)

                % cutoff for FSI/RF 
lowercutoff(1).c = cat(2,0:0.0005:0.2,0.205:0.005:0.5); % for scoringopt = 1
%lowercutoff(2).c = cat(2,0:0.001:0.1,0.11:0.01:1);     % for scoringopt = 2

                % signal to noise ratio
sn(1).c = 2;    % for scoringopt = 1 (cat(2,0:0.2:4,5:1:19,20:5:45,50:10:130); cat(2,0:0.5:19.5,20:5:45,50:10:130);
%sn(2).c = 0;   % for scoringopt = 2

zscores = 0;    % 0: z-score transformation off
                % 1: z-score transformation on
                
                % database to compare to
database = {'e','eb'};  % e: EcoCyc database
                        % eb: BRENDA and ecocyc database 

% folder for results
resfolder = 'results_ROCcurves';

% names to generate files
fnames = {'FSI','RF'};
fnames_2 = {'scoreopt1','scoreopt2'};
fnames_3 = {'EcoCyc','EcoCyc & BRENDA'};

if ~exist(resfolder, 'dir')
    mkdir(resfolder); % folder for output
    disp('Created folder for ROC output');
end
%%
   
% perform parameter scan

out = struct;
ctr = 1;

for j = 1:length(fsamode)
    
    main('.\',fsamode(j),0);
    
    disp(char(strcat('fsa=',num2str(fsamode(j)))));
    
    for d = 1:length(database)
        
        %run once to get hitmap
        [~,hitinfo] = Generate_Interaction_Map(0,0,0,database(d),1);
        
        for l = 1:length(scoringopt)
            disp(char(strcat('scoringopt=',num2str(scoringopt(l)))));
            
            for i = 1:length(sn(scoringopt(l)).c)
                disp(char(strcat('sn=',num2str(sn(scoringopt(l)).c(i)))));
                
                for t = 1:length(zscores)
                    
                    for k = 1:length(lowercutoff(scoringopt(l)).c)
                        [stats,~] = Generate_Interaction_Map(zscores(t),lowercutoff(scoringopt(l)).c(k),sn(scoringopt(l)).c(i),database(d),scoringopt(l), hitinfo);
                        close all hidden;
                        % save stats
                        out(ctr).run = ctr;
                        out(ctr).sn = sn(scoringopt(l)).c(i);
                        out(ctr).fsa = fsamode(j);
                        out(ctr).database = database(d);
                        out(ctr).scoringopt = scoringopt(l);
                        out(ctr).zscore = zscores(t);
                        out(ctr).lowercutoff = lowercutoff(scoringopt(l)).c(k);
                           
                        out(ctr).Datapoints = stats.numofvalues;
                        out(ctr).DatapointsSign = stats.valabovezero;
                        out(ctr).DatapointsSignPerc = stats.valabovezeroperc;
                        
                        out(ctr).totalknown = stats.totalknown;
                        out(ctr).knownDetected = stats.detected;
                        out(ctr).knownNotdetected = stats.notdetected;
                        out(ctr).DetectedPercent = stats.detectedperc ;
                        out(ctr).RatioDetectedvsSign = stats.detectedperc/stats.valabovezeroperc ;
                        out(ctr).eval = out(ctr).RatioDetectedvsSign * out(ctr).knownDetected;
                        
                        % for plotting ROC curves
                        out(ctr).TPR = stats.detected/stats.totalknown;
                        out(ctr).FPR = (stats.valabovezero-stats.detected)/(stats.numofvalues-stats.totalknown);
                        
                        ctr = ctr + 1;
                    end
                end
            end
        end
    end
end

% save all results in matlab structure
filename = fullfile(resfolder,'ROC_Statistics.mat');
save(filename,'out');

%% plot ROC curves

bestparams = struct;    % for parameters that are closest to FPR of 0.05
AUC = struct;           % for info about AUCs

% for temp storage
maxAUCall = max([out.TPR])*max([out.FPR]); 
Actr = 1;
ctr = 1;
names = {};
TP = [];
FP = [];

for i = fsamode
   
    for l = scoringopt
        
        for d = 1:length(database)
            
            % extract relevant data
            datasub = out([out.fsa]==i);
            datasub = datasub([datasub.scoringopt]==l);
            datasub = datasub(strcmp([datasub.database],database(d)));
            
            names(end+1) = {strcat(fnames{i},',',fnames_3{d})};
            fname = char(strcat(fnames{i},'_',fnames_2{l},'_',fnames_3{d},'_'));
            
            
            for j = sn(l).c % one for every scanned S/N ratop
                
                % extract relevant data
                datasubs = datasub([datasub.sn]==j);

                TP(end+1,:) = [datasubs.TPR];
                FP(end+1,:) = [datasubs.FPR];
                
                % save info in structure
                AUC(Actr).fsamode = fnames{i};
                AUC(Actr).scoringopt = l;
                AUC(Actr).sn = j;
                if l == 1 % save type of scoring option that was used
                    AUC(Actr).mode = 'S/N';
                else
                    AUC(Actr).mode = 'combined';
                end
                AUC(Actr).database = fnames_3{d};
                AUC(Actr).AUCraw = abs(trapz(FP(end,:),TP(end,:)));
                AUC(Actr).AUC = abs(trapz(FP(end,:),TP(end,:)))/maxAUCall;
                Actr = Actr +1;
                
                % find parameters that are closest to FPR of 0.05
                bestparams(ctr).fsa = fnames{i};
                bestparams(ctr).scoringopt = l;
                bestparams(ctr).sn = j;
                
                allpos = [datasubs.FPR]-0.05;
                allpos(allpos>0) = 1;
                [~,bestpos] = min(abs(allpos));
                
                bestparams(ctr).FPR = datasubs(bestpos).FPR;
                bestparams(ctr).TPR = datasubs(bestpos).TPR;
                bestparams(ctr).lowercutoff = datasubs(bestpos).lowercutoff;
                
                ctr = ctr+1;
                
            end
        end
    end
end

% save ROC curve 
filename = fullfile(resfolder,'ROC.pdf');
titlename = 'ROC curve';
plotROC(filename,names,TP,FP,titlename);
        
% save best params
filename = fullfile(resfolder,'Best_parameters.mat');
save(filename,'bestparams');

% save AUC info
filename = fullfile(resfolder,'AUC.mat');
save(filename,'AUC');

% save AUC pic
filename = fullfile(resfolder,'AUC.png');
plotAUC(AUC,filename);

end


function plotROC(filename,names,TP,FP,titlename)

f = figure('visible','off');
hold all;
set(gca,'FontSize',6);
colormatrix = colormap(jet(length(names)));


% figure out highest value
maxval = 0;
for i = 1:length(names)
   maxval = max([maxval FP(i,:) TP(i,:)]);
end

for i = 1:length(names)
    color = colormatrix(i,:);
    p(i) = plot(FP(i,:),TP(i,:),'Color',color,'LineWidth',1);
end

% plot line on diagonal
plot([0 1],[0 1],'k--');

ylabel('TPR');
xlabel('FPR');
title(titlename);


ylim([0 maxval]);
xlim([0 maxval]);
legend([p],[names],'location','southeast');

hold off;
saveas(f,filename);

end


function plotAUC(AUC,filename)

f = figure('visible','off');
hold all;
set(gca,'FontSize',8);

bar([AUC.AUC]);

names = {};
for k = 1:length(AUC)
    names{k} = strcat(AUC(k).fsamode,',S/N:',num2str(AUC(k).sn),',db:',AUC(k).database);
end

text(1:length(AUC),[AUC.AUC],num2str([AUC.AUC]'),'vert','bottom','horiz','center'); 
box off;

xticks(1:length(names));
xticklabels(names);
xtickangle(90);

ylabel('tested setups');
xlabel('AUC');
hold off;
saveas(f,filename);

end


