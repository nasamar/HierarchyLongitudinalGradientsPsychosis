%% Script to plot brain maps of cortical hierarchy and SCZ epicenters, and their
%% correlation with MIND associations.

% Copyright (C) 2026 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Hierarchy Longitudinal Gradients Psychosis toolkit.
%
% Hierarchy Longitudinal Gradients Psychosis toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Hierarchy Longitudinal Gradients Psychosis toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Hierarchy Longitudinal Gradients Psychosis toolkit. If not, see 
% <https://www.gnu.org/licenses/>.

close all
clear

location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\';

normalization = '\CN';
% normalization = '\FEP\CN';
normalization = '\FEP+CN';

residual_or_raw = 'raw';
% residual_or_raw = 'residuals';

parcellation = 'aparc_500_sym';
fsa = 'fsa';

degree_68_FEP = readtable([location,'Code\3. MIND_long\Data\degree\',parcellation,'\COMBATLS_covars\degree_68_FEP_residuals.csv'],ReadRowNames=true);
degree_68_CN = readtable([location,'Code\3. MIND_long\Data\degree\',parcellation,'\COMBATLS_covars\degree_68_CN_residuals.csv'],ReadRowNames=true);

gradients_G1_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\residuals\gradients_G1_FEP.csv'],"ReadRowNames",true);
gradients_G2_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\residuals\gradients_G2_FEP.csv'],"ReadRowNames",true);
gradients_G1_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\residuals\gradients_G1_CN.csv'],"ReadRowNames",true);
gradients_G2_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\residuals\gradients_G2_CN.csv'],"ReadRowNames",true);

molecular_maps = readtable([location,'Molecular\Predictive_Maps\all_microsc_DesikanKilliany68.csv'],ReadVariableNames=true);

volumes = readtable([location,'Code\2. MIND\Data\volumes\CorticalMeasuresENIGMA_GrayAvg.csv'],'ReadRowNames',true);

% ENIGMA epicenters
if ~exist([location,'datasets\ENIGMA\fc_ctx_epi.mat'])
    sum_stats = load_summary_stats('schizophrenia');
    
    % Load cortico-cortical functional connectivity data
    [fc_ctx, fc_ctx_labels, ~, ~] = load_fc();
    
    % Load cortico-cortical structural connectivity data
    [sc_ctx, sc_ctx_labels, ~, ~] = load_sc();
    
    % Plot cortico-cortical connectivity matrices
    % f = figure,
    %   imagesc(fc_ctx, [0 0.8]);
    %   axis square;
    %   colormap(Reds);
    %   colorbar;
    %   set(gca, 'YTick', 1:1:length(fc_ctx_labels), ...
    %       'YTickLabel', fc_ctx_labels);
    % 
    % f = figure,
    %   imagesc(sc_ctx, [0 10]);
    %   axis square;
    %   colormap(Blues);
    %   colorbar;
    %   set(gca, 'YTick', 1:1:length(sc_ctx_labels), ...
    %       'YTickLabel', sc_ctx_labels);
    
    CT = sum_stats.CortThick_case_vs_controls;
    SA = sum_stats.CortSurf_case_vs_controls;
    
    % Extract Cohen's d values
    CT_d_DK = CT.d_icv;
    SA_d = SA.d_icv;
    
    
    % % META-ANALYSIS
    % Identify cortical epicenter values (from functional connectivity)
    fc_ctx_epi_DK = zeros(size(fc_ctx, 1), 1);
    fc_ctx_epi_p = zeros(size(fc_ctx, 1), 1);
    for seed = 1:size(fc_ctx, 1)
        seed_conn = fc_ctx(:, seed);
        r_tmp = corrcoef(seed_conn, CT_d_DK);
        fc_ctx_epi_DK(seed) = r_tmp(1, 2);
        fc_ctx_epi_p(seed) = spin_test(seed_conn, CT_d_DK, 'surface_name', 'fsa5', 'parcellation_name', ...
            'aparc', 'n_rot', 1000, 'type', 'pearson');
    end
    
    % Identify cortical epicenter values (from structural connectivity)
    sc_ctx_epi_DK = zeros(size(sc_ctx, 1), 1);
    sc_ctx_epi_p = zeros(size(sc_ctx, 1), 1);
    for seed = 1:size(sc_ctx, 1)
        seed_conn = sc_ctx(:, seed);
        r_tmp = corrcoef(seed_conn, CT_d_DK);
        sc_ctx_epi_DK(seed) = r_tmp(1, 2);
        sc_ctx_epi_p(seed) = spin_test(seed_conn, CT_d_DK, 'surface_name', 'fsa5', 'parcellation_name', ...
            'aparc', 'n_rot', 1000, 'type', 'pearson');
    end
    
    % % Identify cortical epicenter values (from functional connectivity)
    % fc_sax_epi = zeros(size(fc_ctx, 1), 1);
    % fc_sax_epi_p = zeros(size(fc_ctx, 1), 1);
    % for seed = 1:size(fc_ctx, 1)
    %     seed_conn = fc_ctx(:, seed);
    %     r_tmp = corrcoef(seed_conn, SA_d);
    %     fc_sax_epi(seed) = r_tmp(1, 2);
    %     fc_sax_epi_p(seed) = spin_test(seed_conn, SA_d, 'surface_name', 'fsa5', 'parcellation_name', ...
    %         'aparc', 'n_rot', 1000, 'type', 'pearson');
    % end
    % 
    % % % Identify cortical epicenter values (from structural connectivity)
    % sc_sax_epi = zeros(size(sc_ctx, 1), 1);
    % sc_sax_epi_p = zeros(size(sc_ctx, 1), 1);
    % for seed = 1:size(sc_ctx, 1)
    %     seed_conn = sc_ctx(:, seed);
    %     r_tmp = corrcoef(seed_conn, SA_d);
    %     sc_sax_epi(seed) = r_tmp(1, 2);
    %     sc_sax_epi_p(seed) = spin_test(seed_conn, SA_d, 'surface_name', 'fsa5', 'parcellation_name', ...
    %         'aparc', 'n_rot', 1000, 'type', 'pearson');
    % end
    % 
    % Project the results on the surface brain
    % Selecting only regions with p < 0.1 (functional epicenters)
    % fc_ctx_epi_p_sig = zeros(length(fc_ctx_epi_p), 1);
    % fc_ctx_epi_p_sig(find(fc_ctx_epi_p < 0.05)) = fc_ctx_epi_DK(fc_ctx_epi_p<0.05);
    % f = figure,
    %     plot_cortical(parcel_to_surface(fc_ctx_epi_p_sig, 'aparc_fsa5'), ...
    %                 'color_range', [-0.5 0.5], 'cmap', 'GyRd_r')
    % 
    % % Selecting only regions with p < 0.1 (structural epicenters)
    % sc_ctx_epi_p_sig = zeros(length(sc_ctx_epi_p), 1);
    % sc_ctx_epi_p_sig(find(sc_ctx_epi_p < 0.05)) = sc_ctx_epi_DK(sc_ctx_epi_p<0.05);
    % f = figure,
    %     plot_cortical(parcel_to_surface(sc_ctx_epi_p_sig, 'aparc_fsa5'), ...
    %                 'color_range', [-0.5 0.5], 'cmap', 'GyBu_r')
    
    
    save([location,'datasets\ENIGMA\fc_ctx_epi.mat'],"fc_ctx_epi_DK")
    save([location,'datasets\ENIGMA\fc_ctx_epi_p.mat'],"fc_ctx_epi_p")
    save([location,'datasets\ENIGMA\sc_ctx_epi.mat'],"sc_ctx_epi_DK")
    save([location,'datasets\ENIGMA\sc_ctx_epi_p.mat'],"sc_ctx_epi_p")
    
else
    load([location,'datasets\ENIGMA\fc_ctx_epi.mat'])
    load([location,'datasets\ENIGMA\sc_ctx_epi.mat'])

end


[~,idx_to_alphabeth] = sort(volumes(:,1:68).Properties.VariableNames);


aparc_TO_500_sym_names = regexprep(degree_68_FEP.Properties.VariableNames, '^([^_]*_[^_]*)_.*', '$1');
[~,idx] = ismember(aparc_TO_500_sym_names,unique(aparc_TO_500_sym_names));


fc_ctx_epi_DK_ordered = fc_ctx_epi_DK(idx_to_alphabeth);
sc_ctx_epi_DK_ordered = sc_ctx_epi_DK(idx_to_alphabeth);
for i = 1:height(sc_ctx_epi_DK_ordered)
    regions = find(idx==i);
    
    fc_ctx_epi(regions(1):regions(1)+length(regions)-1,:) = fc_ctx_epi_DK_ordered(i,1);

    sc_ctx_epi(regions(1):regions(1)+length(regions)-1,:) = sc_ctx_epi_DK_ordered(i,1);

end

names_DK = readtable([location,'Molecular\parcellations\aparc_500_sym\500.sym_names.txt'],"ReadVariableNames",false);
names_DK = cellfun(@(x, y, z) strcat(x,'_',y,'_', z), names_DK{:,1},names_DK{:,2},names_DK{:,3}, 'UniformOutput', false);    
[~,idx_to_ENIGMA] = ismember(names_DK,degree_68_FEP.Properties.VariableNames);
[~,idx_to_alphabeth] = ismember(degree_68_FEP.Properties.VariableNames,names_DK);


X_rotated = readtable('C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Molecular\parcellations\rotation\500.sym_rotated.txt');
X_rotated = table2array(X_rotated(17:end,:));
X_rotated = X_rotated - 16;


fc_ctx_epi = fc_ctx_epi(idx_to_ENIGMA);
sc_ctx_epi = sc_ctx_epi(idx_to_ENIGMA);


dir_maps = 'C:\Users\usuario\neuromaps-data\annotations';
folders = dir(dir_maps);
folders = folders(3:end);
neuromaps = table();
n_maps = 1;
for k = 1:length(folders)
    subdir_maps = fullfile(dir_maps, folders(k).name);
    
    sub_folders_maps = dir(subdir_maps);
    sub_folders_maps = sub_folders_maps(3:end);

    for l = 1:length(sub_folders_maps)
        % Verify whether 'fsaverage' folder exists    
        fsavg_path = fullfile(fullfile(subdir_maps, sub_folders_maps(l).name), 'fsaverage');
        if exist(fsavg_path, 'dir') 

            folders_fsaverage = dir(fsavg_path);
            folders_fsaverage = folders_fsaverage(3:end);
            
            maps = [];
            for hemisphere = 1:length(folders_fsaverage)
                map_gii = gifti(fullfile(fsavg_path,folders_fsaverage(hemisphere).name));
                maps = [maps;map_gii.cdata];

                if length(folders_fsaverage) < 2
                    maps = [map_gii.cdata;map_gii.cdata];
                end
                
            end

            map = surface_to_parcel(maps',[parcellation,'_',fsa]);
            map = zscore(map(2:end));

            if length(folders_fsaverage) < 2
                map = [map(length(map)/2+1:end),map(length(map)/2+1:end)];
            end

            name_map = sub_folders_maps(l).name;
            
            if contains(name_map,{'SAaxis','evoexp','fcgradient01'})
                % figure('Position', [488   242   560  200])
                plot_cortical( ...
                    parcel_to_surface(map, [parcellation,'_',fsa]), ...
                    'surface_name', fsa, ...
                    'color_range',[-2.97 1.88], 'label_text',name_map,...
                    'position_colorbar','East')
                    % 'color_range',[-max(abs(map)) max(abs(map))],'label_text',name_map)
                colorbar_white_centered([-2.97 1.88],0,7,6)
            end

            neuromaps(n_maps,:) = array2table(map,'VariableNames',names_DK);
            neuromaps.Properties.RowNames(n_maps) = {name_map}; 

        end
        n_maps = n_maps + 1;
    end
end

neuromaps = [neuromaps;array2table([fc_ctx_epi,sc_ctx_epi]',"RowNames",{'fc ctx epi','sc ctx epi'},"VariableNames",neuromaps.Properties.VariableNames)];

for var = {'fc ctx epi','sc ctx epi'}
    map = neuromaps{var,:};
    plot_cortical( ...
        parcel_to_surface(map, [parcellation,'_',fsa]), ...
        'surface_name', fsa, ...
        'color_range',[-0.64 0.26], 'label_text',var{:}, ...
        'position_colorbar','East')

        % 'color_range',[min(map) max(map)], 'label_text',var{:}
        % 'position_colorbar','East')
    colorbar_white_centered([-0.64 0.26],0,7,6)
end


X_rot_ordered = X_rotated(idx_to_alphabeth,:)';


% Interpolation results
nspin = 10000;
r_matrix_CN_baseline = [];
p_spin_matrix_CN_baseline = [];
r_partial_matrix_baseline = [];
p_spin_matrix_baseline = [];
% for gradient = {'','_G2'}
for gradient = {'','_G1','_G2'}
    if strcmp(gradient,'')
        var = 'degrees';
    else
        var = 'gradients';
    end
    for period = {''}
    % for period = {'_baseline',''}
        for cognition = {''}
        % for cognition = {'','_BPRS'}

            variable = [gradient{:},period{:},cognition{:}];
            
            if ~contains(variable,{'BPRS'}) 
                if contains(variable,'baseline')
                    x_vars = {'dx1'};
                    x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler'};
                    
                else
                    x_vars = {'dx1','Treatment_Time','dx1_Treatment_Time','CPZ_equivalent','CPZ_equivalent_Treatment_Time'};
                    % x_vars = {'Treatment_Time','CPZ_equivalent','CPZ_equivalent_Treatment_Time'}; % FEP
                    % x_vars = {'Treatment_Time'}; % CN
                    x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler','protocol_Change15'};

                end
            
            elseif contains(variable,'baseline') 
                x_vars = {['global_',var]};

            else
                x_vars = {['',var],'Treatment_Time',[var,'_Treatment_Time'],'CPZ_equivalent','CPZ_equivalent_Treatment_Time'};

            end

            if strcmp(var,'degrees')
                opts = detectImportOptions([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],'ReadVariableNames',true);
                opts = setvartype(opts, strcat('res_',x_vars), 'double'); 
                Interpolation_Results = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],opts);

            else
                opts = detectImportOptions([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_',var,variable,'.csv'],'ReadVariableNames',true);
                % opts = detectImportOptions([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_',var,variable,'_FEP.csv'],'ReadVariableNames',true);
                opts = setvartype(opts, strcat('res_',x_vars), 'double'); 
                Interpolation_Results = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_',var,variable,'.csv'],opts);
                % Interpolation_Results = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_',var,variable,'_FEP.csv'],opts);
            end

            % poolobj = gcp('nocreate');
            % delete(poolobj);
            % parpool(length(x_vars))


            if contains(variable,'_G1')

                var_partial = mean(gradients_G1_CN{contains(gradients_G1_CN.Properties.RowNames,'_001'),:})';
                var_name = 'G1 CN';

                % var_partial = mean(gradients_G1_FEP{contains(gradients_G1_FEP.Properties.RowNames,'_001'),:})';
                % var_name = 'G1 FEP';

            elseif contains(variable,'_G2')
                
                var_partial = mean(gradients_G2_CN{contains(gradients_G2_CN.Properties.RowNames,'_001'),:})';
                var_name = 'G2 CN';

                % var_partial = mean(gradients_G2_FEP{contains(gradients_G2_FEP.Properties.RowNames,'_001'),:})';
                % var_name = 'G2 FEP';

            else
                var_partial = mean(degree_68_CN{contains(degree_68_CN.Properties.RowNames,'_001'),:})';
                var_name = 'degree CN';

                % var_partial = mean(degree_68_FEP{contains(degree_68_FEP.Properties.RowNames,'_001'),:})';
                % var_name = 'degree FEP';
                
            end

            lim_max = max(max(Interpolation_Results{:,cellfun(@(x) ['res_complete_' x], x_vars, 'UniformOutput', false)}));
            lim_min = min(min(Interpolation_Results{:,cellfun(@(x) ['res_complete_' x], x_vars, 'UniformOutput', false)}));
            % lim_max = 10.13;
            % lim_min = -8.6;

            if strcmp(period,'_baseline') 
                if ~contains(cognition,{'_BPRS'})
                    lim_max = 4;
                    lim_min = -4.1;
                else
                    lim_max = 5.14;
                    lim_min = -6.06;
                end
            end

            r_partial_matrix_tot = table;
            p_partial_matrix_tot = table;
            r_matrix_tot = table;
            p_spin_matrix_tot = table;
            r_matrix_CN_tot = table;
            p_spin_matrix_CN_tot = table;
            for i = 1:length(x_vars)
    
                x_var = x_vars(i);


                % Interpolation_Results{isnan(Interpolation_Results{:,['res_',x_var{:}]}),['res_',x_var{:}]} = 0;
                
                % if contains(x_var{:},'dx1') && ~contains(x_var{:},'Treatment_Time')
                % if ~contains(x_var{:},'CPZ_equivalent_Treatment_Time')
                % if contains(x_var{:},'CPZ') ||  contains(x_var{:},'Treatment_Time')
                % if true
                    
                    if strcmp(period,'_baseline')
                    % if strcmp(x_var,x_vars(1))
                        r = corr(Interpolation_Results{:,['res_complete_',x_var{:}]},var_partial);
           
                        for ispin = 1:nspin
                            betha_perm = Interpolation_Results{X_rot_ordered(ispin,:),['res_complete_',x_var{:}]};
                            r_spin(ispin) = corr(betha_perm,var_partial);
                        end
                        pval_spin = sum(abs(r_spin) > abs(r))/nspin;
                        
    
                        figure;
                        scatter(Interpolation_Results{:,['res_complete_',x_var{:}]},var_partial,'filled')
                        h = lsline;
                        set(h, 'Color', [0 0.4470 0.7410], 'LineWidth', 2);
                        xlabel(['β-value ',x_var{:,:}])
                        ylabel(var_name)
                        title(replace(variable,'_',' '))
                        set(gca, 'FontSize', 18)
                        ylimVals = ylim;
                        yticks(linspace(ylimVals(1), ylimVals(2), 4)); % mostrar al menos 4 nums
                        ytickformat('%.2f');
                        if pval_spin < 0.05
                            annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; ','\color{red}','P_{spin} = ',num2str(pval_spin,2)],'EdgeColor', 'none', 'FontSize', 12);
                        else
                            annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; ','P_{spin} = ',num2str(pval_spin,2)],'EdgeColor', 'none', 'FontSize', 12);
                        end    
                    end

                    r_matrix = table();
                    r_partial_matrix = table();
                    p_partial_matrix = table();
                    p_spin_matrix = table();
                    r_matrix_CN = table();
                    p_spin_matrix_CN = table();
                    for neuromap = 1:height(neuromaps)

                        % if contains(neuromaps.Properties.RowNames(neuromap),{'SAaxis','fcgradient01','evoexp','myelin'}) || corr(Interpolation_Results{:,['res_complete_',x_var{:}]},neuromaps{neuromap,idx}') >= 0.3
                        if contains(neuromaps.Properties.RowNames(neuromap),{'SAaxis','fcgradient01','evoexp','fc ctx epi','sc ctx epi'})
                        % if ~contains(neuromaps.Properties.RowNames(neuromap),{'SAaxis','fcgradient01','evoexp','myelin'}) && partialcorr(Interpolation_Results{:,['res_complete_',x_var{:}]},neuromaps{neuromap,idx}',var_partial) >= 0.3

                           [r_partial,p_partial] = partialcorr(Interpolation_Results{:,['res_complete_',x_var{:}]},neuromaps{neuromap,idx_to_alphabeth}',var_partial);
                            r = corr(Interpolation_Results{:,['res_complete_',x_var{:}]},neuromaps{neuromap,idx_to_alphabeth}');

                            for ispin = 1:nspin
                                betha_perm = Interpolation_Results{X_rot_ordered(ispin,:),['res_complete_',x_var{:}]};
                                r_spin(ispin) = partialcorr(betha_perm,neuromaps{neuromap,idx_to_alphabeth}',var_partial);
                            end
                            

                            pval_spin = sum(abs(r_spin) > abs(r_partial))/nspin;

                            
                            % figure;
                            % scatter(Interpolation_Results{:,['res_complete_',x_var{:}]},neuromaps{neuromap,idx_to_alphabeth},40,var_partial,'filled')
                            % lsline
                            % xlabel(['β-value ',x_var{:,:}])
                            % ylabel(neuromaps.Properties.RowNames(neuromap))
                            % if pval_spin < 0.05
                            %     annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; r_{partial} = ',num2str(r_partial,2),'; ','\color{red}','P_{spin} = ',num2str(pval_spin,2),'\color{black}'],'EdgeColor', 'none', 'FontSize', 12);
                            % else
                            %     annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; r_{partial} = ',num2str(r_partial,2),'; P_{spin} = ',num2str(pval_spin,2),'\color{red}'],'EdgeColor', 'none', 'FontSize', 12);
                            % end
                            % title(replace(variable,'_',' '))

                            r_matrix = [r_matrix;array2table(r,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",x_var)];
                            r_partial_matrix = [r_partial_matrix;array2table(r_partial,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",x_var)];
                            p_partial_matrix = [p_partial_matrix;array2table(p_partial,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",x_var)];
                            p_spin_matrix = [p_spin_matrix;array2table(pval_spin,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",x_var)];

                            if strcmp(x_var,x_vars(1))
                            % if strcmp(x_var,x_vars(1)) & strcmp(period,'_baseline')

                                r = corr(var_partial,neuromaps{neuromap,idx_to_alphabeth}');
                                for ispin = 1:nspin
                                    var_partial_perm = var_partial(X_rot_ordered(ispin,:));
                                    r_spin(ispin) = corr(var_partial_perm,neuromaps{neuromap,idx_to_alphabeth}');
                                end
                                pval_spin = sum(abs(r_spin) > abs(r))/nspin;
    
                                
                                % figure;
                                % scatter(var_partial,neuromaps{neuromap,idx_to_alphabeth},'filled')
                                % lsline
                                % xlabel(var_name)
                                % ylabel(neuromaps.Properties.RowNames(neuromap))
                                % if pval_spin < 0.05
                                %     annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; ','\color{red}','P_{spin} = ',num2str(pval_spin,2)],'EdgeColor', 'none', 'FontSize', 12);
                                % else
                                %     annotation('textbox', [0.15, 0.25, 0.3, 0.05], 'String', ['r = ',num2str(r,2),'; P_{spin} = ',num2str(pval_spin,2)],'EdgeColor', 'none', 'FontSize', 12);
                                % end
                                % title(replace(variable,'_',' '))

                                r_matrix_CN = [r_matrix_CN;array2table(r,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",{var_name})];
                                p_spin_matrix_CN = [p_spin_matrix_CN;array2table(pval_spin,"RowNames",neuromaps.Properties.RowNames(neuromap),"VariableNames",{var_name})];

                            end
                        end
                    end 

                    r_matrix_tot = [r_matrix_tot,r_matrix];
                    r_partial_matrix_tot = [r_partial_matrix_tot,r_partial_matrix];
                    p_partial_matrix_tot = [p_partial_matrix_tot,p_partial_matrix];
                    p_spin_matrix_tot = [p_spin_matrix_tot,p_spin_matrix];
                    r_matrix_CN_tot = [r_matrix_CN_tot,r_matrix_CN];
                    p_spin_matrix_CN_tot = [p_spin_matrix_CN_tot,p_spin_matrix_CN];
                    
                % end
            end

            % if false
            if strcmp(period,'')
                % figure;
                % heatmap(r_matrix_tot{:,:},'XData',r_matrix_tot.Properties.VariableNames,'YData',r_matrix_tot.Properties.RowNames,'Title',['r',replace(variable,'_',' ')]);
                % colorbar_white_centered(r_matrix_tot{:,:})
    
                figure;
                heatmap(r_partial_matrix_tot{:,:},'XData',r_partial_matrix_tot.Properties.VariableNames,'YData',r_partial_matrix_tot.Properties.RowNames,'Title',['r partial',replace(variable,'_',' ')],CellLabelFormat='%.2f');
                colorbar_white_centered(r_partial_matrix_tot{:,:})
                
                figure;
                heatmap(r_partial_matrix_tot{:,:}','YData',r_partial_matrix_tot.Properties.VariableNames,'XData',r_partial_matrix_tot.Properties.RowNames,'Title',['r partial',replace(variable,'_',' ')],CellLabelFormat='%.2f');
                colorbar_white_centered(r_partial_matrix_tot{:,:})

                % figure;
                % heatmap(p_partial_matrix_tot{:,:},'XData',r_partial_matrix_tot.Properties.VariableNames,'YData',r_partial_matrix_tot.Properties.RowNames,'Title',['p partial',replace(variable,'_',' ')]);
                % colorbar_white_centered(p_partial_matrix_tot{:,:})
    
                % figure;
                % heatmap(p_spin_matrix{:,:},'XData',r_matrix.Properties.VariableNames,'YData',r_matrix.Properties.RowNames,'Title',['p spin',replace(variable,'_',' ')]);
    
                for j = 1:width(p_spin_matrix_tot{:,:})
                    p_spin_matrix_tot_corrected(:,j) = mafdr(p_spin_matrix_tot{:,j},'BHFDR','true');
                end
                
                figure;
                heatmap(p_spin_matrix_tot_corrected,'XData',r_partial_matrix_tot.Properties.VariableNames,'YData',r_partial_matrix_tot.Properties.RowNames,'Title',['p spin',replace(variable,'_',' ')]);
                
            end
        end

        if strcmp(period,'_baseline')
            r_matrix_CN_baseline = [r_matrix_CN_baseline,r_matrix_CN{:,:}];
            p_spin_matrix_CN_baseline = [p_spin_matrix_CN_baseline,p_spin_matrix_CN{:,:}];
            r_partial_matrix_baseline = [r_partial_matrix_baseline,r_partial_matrix{:,:}];
            p_spin_matrix_baseline = [p_spin_matrix_baseline,p_spin_matrix{:,:}];
        end

    end
end

figure;
h=heatmap(r_partial_matrix_tot{:,:}','YData',r_partial_matrix_tot.Properties.VariableNames,'XData',r_partial_matrix_tot.Properties.RowNames,'Title',['r partial',replace(variable,'_',' ')],CellLabelFormat='%.2f');
colorbar_white_centered(r_partial_matrix_tot{:,:})
clim = h.ColorLimits;
ax = axes('Position',[0.05 -0.002 0.90 1],'Visible','off','FontSize',20,'TickDir','out');
colorbar(ax,'southoutside','Orientation','horizontal');
caxis(ax, clim);

figure;
h = heatmap(r_matrix_CN_baseline,'XData',{'degrees','G1','G2'},'YData',r_partial_matrix.Properties.RowNames,'Title','r','FontSize',20,CellLabelFormat='%.2f');
colorbar_white_centered(r_matrix_CN_baseline)
% Horizontal colorbar
clim = h.ColorLimits;
ax = axes('Position',[0.05 -0.002 0.90 1],'Visible','off','FontSize',20,'TickDir','out');
colorbar(ax,'southoutside','Orientation','horizontal');
caxis(ax, clim);


for i = 1:size(p_spin_matrix_CN_baseline,2)
    p_spin_matrix_CN_baseline_corrected(:,i) = mafdr(p_spin_matrix_CN_baseline(:,i),'BHFDR','true');
end

figure;
heatmap(p_spin_matrix_CN_baseline_corrected','XData',r_partial_matrix.Properties.RowNames,'YData',{'degrees','G1','G2'},'Title','p spin');


for i = 1:size(p_spin_matrix_baseline,2)
    p_spin_matrix_baseline_corrected(:,i) = mafdr(p_spin_matrix_baseline(:,i),'BHFDR','true');    
end

figure;
heatmap(r_partial_matrix_baseline','XData',r_partial_matrix.Properties.RowNames,'YData',{'dx1 degrees','dx1 G1','dx1 G2'},'Title','r partial','FontSize',20,CellLabelFormat='%.2f');
colorbar_white_centered(r_partial_matrix_baseline)

figure;
heatmap(p_spin_matrix_baseline_corrected','XData',r_partial_matrix.Properties.RowNames,'YData',{'dx1 degrees','dx1 G1','dx1 G2'},'Title','p spin');



