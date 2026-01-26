%% Script to compute HC and SSD MIND edges and degrees.

% Copyright (C) 2026 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Hierarchy MIND Psychosis toolkit.
%
% Hierarchy MIND Psychosis toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Hierarchy MIND Psychosis toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Hierarchy MIND Psychosis toolkit. If not, see 
% <https://www.gnu.org/licenses/>.

clear
close all
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\';
% type = 'COMBAT';
% type = 'COMBATLS';
type = 'COMBATLS_covars';

% type = '';

parcellation = 'aparc_500_sym';
parcellation = 'subcortical';

edge_metrics = {'edge'};
regional_metrics = {'degree'};

if ~strcmp(parcellation,'aparc')
    edge_metrics = {'edge'};
    regional_metrics = {'degree'};
end

metrics = [edge_metrics,regional_metrics];

for i = 1:length(metrics)
    writing_path = ['Code\3. MIND_long\Data\',metrics{i},'\',parcellation,'\',type,'\'];
    eval([['write_',metrics{i}] ' = writing_path;']);
end

write_edge = ['Code\3. MIND_long\Data\edges\',parcellation,'\',type,'\'];

if contains(type,'COMBATLS')
    read = ['Code\R\Connectivity\',parcellation,'_longitudinal\',type,'\'];
elseif contains(type,'COMBAT')
    read = ['Code\Python\Connectivity\Longitudinal\edges\',parcellation,'\',type,'\'];
else
    read = ['datasets\MIND\MIND_networks_PAFIP\',parcellation,'_longitudinal\'];
end

cd([location,write_degree])


if exist([location,write_degree,'degree_68_CN.csv'])
    for i = 1:length(regional_metrics)
        reading_path_CN = readtable([location,eval(['write_',regional_metrics{i}]),[regional_metrics{i},'_68_CN.csv']],"ReadRowNames",true);
        eval([[regional_metrics{i},'_68_CN'] ' = reading_path_CN;']);
        reading_path_FEP = readtable([location,eval(['write_',regional_metrics{i}]),[regional_metrics{i},'_68_FEP.csv']],"ReadRowNames",true);
        eval([[regional_metrics{i},'_68_FEP'] ' = reading_path_FEP;']);
    end
    for i = 1:length(edge_metrics)
        load([location,eval(['write_',edge_metrics{i}]),[edge_metrics{i},'_68_CN.mat']])
        load([location,eval(['write_',edge_metrics{i}]),[edge_metrics{i},'_68_FEP.mat']])
    end

    
else
    files = dir([location,read]);
    
    files = files(3:end); % all sessions
    
    i_FEP = 1;
    i_CN = 1;
    FEP_subjects = {};
    CN_subjects = {};
    for i = 1:length(files)
        MIND = readtable([[location,read],files(i).name],'ReadRowNames',true,'ReadVariableNames',true);

        % Sort alphabetically
        MIND = sortrows(MIND,'RowNames');
        col_order = sort(MIND.Properties.VariableNames);
        MIND = MIND(:,col_order);
        
        % Averaged hemispheres
        lh = MIND(contains(MIND.Properties.RowNames,'lh_'),contains(MIND.Properties.VariableNames,'lh_'));
        rh = MIND(contains(MIND.Properties.RowNames,'rh_'),contains(MIND.Properties.VariableNames,'rh_'));
        
        try
            mean_MIND = (lh{:,:}+rh{:,:})/2;
        catch
            lacking_region = setdiff(replace(rh.Properties.VariableNames,'rh_',''),replace(lh.Properties.VariableNames,'lh_',''));
            index_region = find(strcmp(rh.Properties.VariableNames,['rh_',lacking_region{:}]));
            lacking_region = ['lh_',lacking_region{:}];
            lh = [lh(:,1:index_region-1), array2table(nan(length(mean_MIND)-1,1),"VariableNames",{lacking_region}), lh(:,index_region:end)];
            lh = [lh(1:index_region-1,:); array2table(nan(1,length(mean_MIND)),"RowNames",{lacking_region},"VariableNames",lh.Properties.VariableNames); lh(index_region:end,:)];
            lh(:,index_region) = rh(:,index_region);
            lh(index_region,:) = rh(index_region,:);       
            mean_MIND = (lh{:,:}+rh{:,:})/2;

            if contains(lacking_region,'lh_')
                MIND = [MIND(:,1:index_region-1), array2table(nan(width(MIND),1),"VariableNames",{lacking_region}), MIND(:,index_region:end)];
                MIND = [MIND(1:index_region-1,:); array2table(nan(1,width(MIND)),"RowNames",{lacking_region},"VariableNames",MIND.Properties.VariableNames); MIND(index_region:end,:)];
                MIND(:,index_region) = MIND(:,width(mean_MIND)+index_region);
                MIND(index_region,:) = MIND(width(mean_MIND)+index_region,:);       
            else
                MIND = [MIND(:,1:width(mean_MIND)+index_region-1), array2table(nan(width(MIND),1),"VariableNames",{lacking_region}), MIND(:,width(mean_MIND)+index_region:end)];
                MIND = [MIND(1:width(mean_MIND)+index_region-1,:); array2table(nan(1,width(MIND)),"RowNames",{lacking_region},"VariableNames",MIND.Properties.VariableNames); MIND(width(mean_MIND)+index_region:end,:)];
                MIND(:,width(mean_MIND)+index_region) = MIND(:,index_region);
                MIND(width(mean_MIND)+index_region,:) = MIND(index_region,:); 
            end

        end
    
        % FEP 
        if str2double(regexp(files(i).name, '(\d+)_', 'tokens', 'once')) < 1000
    
            FEP_subjects = [FEP_subjects,replace(files(i).name,'.csv','')];

            % Edge
            edge_68_FEP{i_FEP,:} = MIND;
    
            % REGIONAL
            % Degree
            degree_68_FEP(i_FEP,:) = mean(MIND{:,:});
            

            i_FEP = i_FEP + 1;
        
        % CN
        else
            CN_subjects = [CN_subjects,replace(files(i).name,'.csv','')];
    
            % Edge
            edge_68_CN{i_CN,:} = MIND;

            % REGIONAL
            % Degree
            degree_68_CN(i_CN,:) = mean(MIND{:,:});
            
            i_CN = i_CN + 1;
        end
    end
    
    if strcmp(type,'')
        CN_subjects = regexprep(strrep(strrep(CN_subjects, 'sub-', ''), 'ses-', ''),'\..*',''); 
        FEP_subjects = regexprep(strrep(strrep(FEP_subjects, 'sub-', ''), 'ses-', ''),'\..*','');
    end

    for edge_metric = edge_metrics
        edge_metric_FEP = [FEP_subjects',eval([edge_metric{:},'_68_FEP'])];
        edge_metric_CN = [CN_subjects',eval([edge_metric{:},'_68_CN'])];
        eval([[edge_metric{:},'_68_FEP'] ' = edge_metric_FEP'])
        eval([[edge_metric{:},'_68_CN'] ' = edge_metric_CN'])

        save([location,eval(['write_',edge_metric{:}]),[edge_metric{:},'_68_FEP.mat']],[edge_metric{:},'_68_FEP'])
        save([location,eval(['write_',edge_metric{:}]),[edge_metric{:},'_68_CN.mat']],[edge_metric{:},'_68_CN'])
    
    end
    for regional_metric = regional_metrics
        regional_metric_FEP = array2table(eval([regional_metric{:},'_68_FEP']),"RowNames",FEP_subjects,"VariableNames",MIND.Properties.VariableNames);
        regional_metric_CN = array2table(eval([regional_metric{:},'_68_CN']),"RowNames",CN_subjects,"VariableNames",MIND.Properties.VariableNames);

        writetable(regional_metric_FEP,[[location,eval(['write_',regional_metric{:}])],['\',regional_metric{:},'_68_FEP.csv']],'WriteRowNames',true);
        writetable(regional_metric_CN,[[location,eval(['write_',regional_metric{:}])],['\',regional_metric{:},'_68_CN.csv']],'WriteRowNames',true);
        save([location,eval(['write_',regional_metric{:}]),[regional_metric{:},'_68_FEP.mat']],[regional_metric{:},'_68_FEP'])
        save([location,eval(['write_',regional_metric{:}]),[regional_metric{:},'_68_CN.mat']],[regional_metric{:},'_68_CN'])
    end

end

degree_long_FEP = degree_68_FEP;
degree_long_CN = degree_68_CN;
degree_long = [degree_68_CN;degree_68_FEP];

edge_long_FEP = edge_68_FEP;
edge_long_CN = edge_68_CN;

edge_long_R_CN = edge_long_CN;
for i = 1:length(edge_68_CN)
    edge_long_R_CN{i,2} = table2array(edge_long_CN{i,2});
end
save([location,write_edge,'edge_long_R_CN.mat'],'edge_long_R_CN')

edge_long_R_FEP = edge_long_FEP;
for i = 1:length(edge_68_FEP)
    edge_long_R_FEP{i,2} = table2array(edge_long_FEP{i,2});
end
save([location,write_edge,'edge_long_R_FEP.mat'],'edge_long_R_FEP')

edge_FEP = edge_long_FEP(cellfun(@(x) ischar(x) && contains(x, '_001'),edge_long_FEP(:,1)),:);
edge_CN = edge_long_CN(cellfun(@(x) ischar(x) && contains(x, '_001'),edge_long_CN(:,1)),:);


% EFFECT SIZES DEGREE
degree_long_FEP.Properties.VariableNames = replace(replace(degree_long_FEP.Properties.VariableNames,'lh','L'),'rh','R');
degree_long_FEP.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_long_FEP.Properties.VariableNames, 'UniformOutput', false);
degree_long_CN.Properties.VariableNames = replace(replace(degree_long_CN.Properties.VariableNames,'lh','L'),'rh','R');
degree_long_CN.Properties.VariableNames = cellfun(@(x) [x, '_grayavg'], degree_long_CN.Properties.VariableNames, 'UniformOutput', false);

etiv = readtable([location,'\datasets\PAFIP\aparc_formatted_1293s_fs_v7.4.1_long_FUP_13-01-2025\COMBATLS\global_aparc_COMBATLS_covars'],"ReadRowNames",true);
etiv = etiv(:,'EstimatedTotalIntraCranialVol');

opts = detectImportOptions([location,'\datasets\PAFIP\Covariates_complete.csv'],"ReadRowNames",true); % Detecta las opciones de importación automáticamente
opts = setvartype(opts, 'Machine_Teslas', 'char');
covariates_long = readtable([location,'\datasets\PAFIP\Covariates_complete.csv'],opts);

covariates_long = covariates_long(etiv.Properties.RowNames,:);

if strcmp(parcellation,'subcortical')
    excluded_subjects = readtable([location,'datasets/PAFIP/ENIGMA_Shape/QA_Status.csv'],ReadRowNames=true);
    excluded_subjects = excluded_subjects(any(excluded_subjects{:,:}==0,2),:);
    excluded_subjects = regexprep(strrep(strrep(excluded_subjects.Properties.RowNames, 'sub-', ''), 'ses-', ''),'\..*',''); 
    covariates_long = covariates_long(~ismember(covariates_long.Properties.RowNames,excluded_subjects),:);
end
degree_long_FEP = degree_long_FEP(ismember(degree_long_FEP.Properties.RowNames,covariates_long.Properties.RowNames),:);
degree_long_CN = degree_long_CN(ismember(degree_long_CN.Properties.RowNames,covariates_long.Properties.RowNames),:);
covariates_long_FEP = covariates_long(degree_long_FEP.Properties.RowNames,:);
covariates_long_CN = covariates_long(degree_long_CN.Properties.RowNames,:);

covariates = [covariates_long_CN;covariates_long_FEP];

covariates_FEP = covariates_long_FEP(covariates_long_FEP.Assessment==1,:);
covariates_CN = covariates_long_CN(covariates_long_CN.Assessment==1,:);

euler_index = readtable([location,'\datasets\PAFIP\euler_longitudinal.csv'],ReadRowNames=true);
euler_index.Properties.RowNames = strrep(strrep(euler_index.Properties.RowNames, 'sub-', ''), 'ses-', ''); 

covariates.euler_lh = euler_index{[degree_long_CN.Properties.RowNames;degree_long_FEP.Properties.RowNames],"euler_lh"};
covariates.euler_rh = euler_index{[degree_long_CN.Properties.RowNames;degree_long_FEP.Properties.RowNames],"euler_rh"};

degree_68_FEP = degree_long_FEP(covariates_long_FEP.Assessment==1,:);
degree_68_CN = degree_long_CN(covariates_long_CN.Assessment==1,:);

etiv_FEP = etiv(degree_long_FEP.Properties.RowNames,:);
etiv_CN = etiv(degree_long_CN.Properties.RowNames,:);


figure;
scatter(mean(degree_long_CN{:,:},2),etiv_CN{:,:},'filled')
corr(mean(degree_long_CN{:,:},2),etiv_CN{:,:})
lsline
xlabel('MIND CN')
ylabel('eTIV')

figure;
scatter(mean(degree_long_CN{:,:},2),euler_index{degree_long_CN.Properties.RowNames,3},'filled')
corr(mean(degree_long_CN{:,:},2),euler_index{degree_long_CN.Properties.RowNames,3})
lsline
xlabel('MIND CN')
ylabel('Euler index')


% Comment for baseline effect (not longitudinal)
% degree_long_FEP = degree_long_FEP(covariates_long_FEP.Assessment==1,:);
% covariates_long_FEP = covariates_long_FEP(covariates_long_FEP.Assessment==1,:);
% etiv_FEP = etiv_FEP(covariates_long_FEP.Assessment==1,:);

% DEGREE
% Site effect
figure;
scatter([1:height(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'GE_1.5T'),:})],mean(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'GE_1.5T'),:},2),'filled')
ylabel(['Global degree ',type])
hold on
scatter([height(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'GE_1.5T'),:})+1:height(degree_long_FEP)],mean(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'Philips_3T'),:},2),'filled')
legend({'PAFIP_1.5T','PAFIP_3T'})
title('FEP')
[h,p] = ttest2(mean(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'GE_1.5T'),:},2),mean(degree_long_FEP{strcmp(covariates_long_FEP.Machine_Teslas,'Philips_3T'),:},2));

figure;
scatter([1:height(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'GE_1.5T'),:})],mean(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'GE_1.5T'),:},2),'filled')
hold on
scatter([height(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'GE_1.5T'),:})+1:height(degree_long_CN)],mean(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'Philips_3T'),:},2),'filled')
legend({'PAFIP_1.5T','PAFIP_3T'})
title('CN')
ylabel(['Global degree ',type])
[h,p] = ttest2(mean(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'GE_1.5T'),:},2),mean(degree_long_CN{strcmp(covariates_long_CN.Machine_Teslas,'Philips_3T'),:},2));


% Sex effect
figure;
scatter([1:height(degree_long_FEP{covariates_long_FEP.Sex==0,:})],mean(degree_long_FEP{covariates_long_FEP.Sex==0,:},2),'filled')
ylabel(['Global degree ',type])
xlabel('SSD patient')
hold on
scatter([height(degree_long_FEP{covariates_long_FEP.Sex==0,:})+1:height(degree_long_FEP)],mean(degree_long_FEP{covariates_long_FEP.Sex==1,:},2),'filled')
legend({'Male','Female'})
title('FEP')
[h,p] = ttest2(mean(degree_long_FEP{covariates_long_FEP.Sex==0,:},2),mean(degree_long_FEP{covariates_long_FEP.Sex==1,:},2));

figure;
scatter([1:height(degree_long_CN{covariates_long_CN.Sex==0,:})],mean(degree_long_CN{covariates_long_CN.Sex==0,:},2),'filled')
hold on
scatter([height(degree_long_CN{covariates_long_CN.Sex==0,:})+1:height(degree_long_CN)],mean(degree_long_CN{covariates_long_CN.Sex==1,:},2),'filled')
legend({'Male','Female'})
title('CN')
ylabel(['Global degree ',type])
[h,p] = ttest2(mean(degree_long_CN{covariates_long_CN.Sex==0,:},2),mean(degree_long_CN{covariates_long_CN.Sex==1,:},2));

figure;
scatter(covariates_long_CN.Age_MRI,mean(degree_long_CN{:,:},2),'filled')
title('CN')
xlabel('Age')
ylabel(['Global degree ',type])
lsline

figure;
scatter(covariates_long_FEP.Age_MRI,mean(degree_long_FEP{:,:},2),'filled')
title('FEP')
xlabel('Age')
ylabel(['Global degree ',type])
lsline

figure;
scatter(etiv_CN{:,:},mean(degree_long_CN{:,:},2),'filled')
title('CN')
xlabel('eTIV')
ylabel(['Global degree ',type])
lsline


figure;
scatter(etiv_FEP{:,:},mean(degree_long_FEP{:,:},2),'filled')
title('FEP')
xlabel('eTIV')
ylabel(['Global degree ',type])
lsline

figure;
scatter([1:height(degree_long_CN{covariates_long_CN.Sex==0,:})],etiv_CN{covariates_long_CN.Sex==0,:},'filled')
hold on
scatter([height(degree_long_CN{covariates_long_CN.Sex==0,:})+1:height(degree_long_CN)],etiv_CN{covariates_long_CN.Sex==1,:},'filled')
legend({'Male','Female'})
title('CN')
ylabel(['eTIV ',type])
[h,p] = ttest2(etiv_CN{covariates_long_CN.Sex==0,:},etiv_CN{covariates_long_CN.Sex==1,:});

figure;
scatter([1:height(degree_long_FEP{covariates_long_FEP.Sex==0,:})],etiv_FEP{covariates_long_FEP.Sex==0,:},'filled')
hold on
scatter([height(degree_long_FEP{covariates_long_FEP.Sex==0,:})+1:height(degree_long_FEP)],etiv_FEP{covariates_long_FEP.Sex==1,:},'filled')
legend({'Male','Female'})
title('FEP')
ylabel(['eTIV ',type])
[h,p] = ttest2(etiv_FEP{covariates_long_FEP.Sex==0,:},etiv_FEP{covariates_long_FEP.Sex==1,:});
