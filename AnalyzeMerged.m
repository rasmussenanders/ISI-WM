%% Plot merged data
clear; clc; close
cd('/Users/med-anr/Google Drive/Documents (GDrive)/Projects/Current/Humanblink short ISI/Analysis')
addpath(genpath('/Users/med-anr/Google Drive/Documents (GDrive)/Projects/Current/Humanblink short ISI/Analysis'))
load('MergedNew.mat') % Load table
%SubjectsNew(42:87,:) = []; % Delete all pilots

%SubjectsNew(,:) = []; % Delete all pilots
%25, 34, 36, 37, 38, 39, 40, 41


%% Loop through subjects
close all
Perc = nan(size(SubjectsNew,1),10);
idx = 1; % Create variable for regression table

for j = 1 : size(SubjectsNew,1)
    
    % Load stuff
    mdmt_sf = SubjectsNew.mdmt_sf{j}; % Collect smoothed data
    codes_cat =  SubjectsNew.codes_cat{j};  % Collect trial codes
    isi = SubjectsNew.ISI(j); % Collect the ISI
    sex = SubjectsNew.Sex(j);
    wm = SubjectsNew.WM(j);
    wmscore = SubjectsNew.WM_Score_Tot(j);
    
    % If there are more than 10 blocks, take the 10 first blocks, if there is less than 10 blocks then take all blocks
    temp =  SubjectsNew.Perc{j};
    if numel(temp) < 10; Perc(j,1:numel(temp)) = temp; else; Perc(j,:) = temp(1:10); temp = temp(1:10); end 
    
    % Expand the table for  regression analysis
    for i = 1: numel(temp) 
        T = table(j,i,sex,wm,wmscore,isi,temp(i));
        
        if idx == 1
            Treg = T;
        else
            Treg = [Treg; T];
        end
        
        idx = idx + 1;
    end
    
end

Treg.Properties.VariableNames = {'FP' 'Block' 'Sex' 'WM' 'WM_Score' 'ISI' 'CRs' };

%% Run One-way random effects model (mixed effects model)
% Run a mixed model (fitlme) to test the effect of the ISI on the
% percentage of CRs. The block number and the ISI will be fixed variables
% since they are determined before the experiment and hence there is no
% sampling error

Treg.Sex = double(Treg.Sex); % Convert to double, Femal = 1; Male = 2
Treg.WM = double(Treg.WM);
%lme = fitlme(Treg, 'CRs ~   ISI + Sex  +  Block + WM  + (1|FP) ')
Sample = find( Treg.ISI == 500);
Treg2 = Treg(Sample,:)
lme = fitlme(Treg, 'CRs ~    Block + Sex + WM + ISI + (1|FP) ')


%% Secon LME with only WM

% TregWM = Treg(~isnan(Treg.WM_Score) & Treg.ISI == 500,:); % Create subtable for rows when there is a value
% lme = fitlme(TregWM, 'CRs ~   ISI + WM_Score + Block + (1|FP) ')

%% Percent CR plot

close all
figure()
ISI150 = Perc(SubjectsNew.ISI==150,:);
ISI250 = Perc(SubjectsNew.ISI==250,:);
ISI500 = Perc(SubjectsNew.ISI==500,:);
    
Sample = ISI150; 
plot(1:10,nanmean(Sample),'ob', 'color', [0 0.4470 0.7410]); hold on
errorbar(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'color', [0 0.4470 0.7410])
xlim([0 11]); 
%boundedline(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'b' ); hold on

Sample = ISI250; 
plot(1:10,nanmean(Sample),'or',  'color',[0.4660 0.6740 0.1880]); hold on
errorbar(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'color', [0.4660 0.6740 0.1880])
xlim([0 11]); 
%boundedline(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'r')




Sample = ISI500; 
plot(1:10,nanmean(Sample),'ok',  'color', [0.6350 0.0780 0.1840]); hold on
errorbar(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'color', [0.6350 0.0780 0.1840])
xlim([0 11]); 
%boundedline(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'k')

ylim([0 1])
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');
h=legend('150',' ' ,'250', '  ' ,'500', 'Location', 'northwest'); set(h,'color','none')

%% WM vs non-WM

figure()

% ISI 500, WM 
Sample = Perc(SubjectsNew.Group == 'Med'  & SubjectsNew.WM == 'Yes'  & SubjectsNew.ISI == 500,:) ;
plot(1:10,nanmean(Sample),'ok',  'color', [0.6350 0.0780 0.1840]); hold on
errorbar(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'color', [0.6350 0.0780 0.1840])
xlim([0 11]); 

% ISI 500 no WM
Sample = Perc(SubjectsNew.Group == 'Med' & SubjectsNew.WM == 'No'  & SubjectsNew.ISI == 500,:) ;
plot(1:10,nanmean(Sample),'ok',  'color', [0.8500 0.3250 0.0980]); hold on
errorbar(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'color', [0.8500 0.3250 0.0980])
xlim([0 11]); 


%title('ISI500, WM (b, n=15), or no WM (r, n=9)')
h=legend('500 WM',' ' ,'500 No WM', '250 WM', '250 no WM', 'Location', 'southeast'); set(h,'color','none')
ylim([0 1])
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');



%% Men vs Women
figure()
Sample = Perc(SubjectsNew.Sex == 'Male'  & SubjectsNew.ISI == 500,:) ;
boundedline(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'b' ); hold on

Sample = Perc(SubjectsNew.Sex == 'Female' & SubjectsNew.ISI == 500,:) ;
boundedline(1:10, nanmean(Sample), nanstd(Sample) / sqrt(size(Sample,1)), 'r')

title('ISI500, Male (b, n=53), Female (r, n=17)')
ylim([0 1])
set(gca, 'TickDir', 'out'); box off; set(gca,'color','none');


