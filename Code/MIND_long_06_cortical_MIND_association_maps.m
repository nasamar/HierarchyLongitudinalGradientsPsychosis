%% Script to plot brain maps of cortical MIND associations.

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

names_DK = readtable([location,'Molecular\parcellations\aparc_500_sym\500.sym_names.txt'],"ReadVariableNames",false);
names_DK = cellfun(@(x, y, z) strcat(x,'_',y,'_', z), names_DK{:,1},names_DK{:,2},names_DK{:,3}, 'UniformOutput', false);    
[~,idx_to_ENIGMA] = ismember(names_DK,degree_68_FEP.Properties.VariableNames);

% Interpolation results
% for gradient = {'','_G2'}
for gradient = {'','_G1','_G2'}
    if strcmp(gradient,'')
        var = 'degrees';
    else
        var = 'gradients';
    end
    % for period = {''}
    for period = {'_baseline',''}
        % for cognition = {''}
        for cognition = {'','_BPRS'}

            variable = [gradient{:},period{:},cognition{:}];
            
            if ~contains(variable,{'BPRS'}) 
                if contains(variable,'baseline')
                    x_vars = {'dx1'};
                    % x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler'};
                    
                else
                    x_vars = {'dx1','Treatment_Time','dx1_Treatment_Time','CPZ_equivalent','CPZ_equivalent_Treatment_Time'};
                    % x_vars = {'Treatment_Time','CPZ_equivalent','CPZ_equivalent_Treatment_Time'}; % FEP
                    % x_vars = {'Treatment_Time'}; % CN
                    % x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler','protocol_Change15'};

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

            
            for i = 1:length(x_vars)
    
                x_var = x_vars(i);


                % Interpolation_Results{isnan(Interpolation_Results{:,['res_',x_var{:}]}),['res_',x_var{:}]} = 0;

                % if contains(x_var{:},'dx1') && ~contains(x_var{:},'Treatment_Time')
                % if contains(variable,'baseline') 
                % if contains(x_var{:},var) 
                    % % figure('Position', [488   242   560  400])
                    % plot_cortical( ...
                    %     [parcel_to_surface(Interpolation_Results{idx_to_ENIGMA,['res_complete_',x_var{:}]}, [parcellation,'_',fsa]);parcel_to_surface(Interpolation_Results{idx_to_ENIGMA,['res_',x_var{:}]}, [parcellation,'_',fsa])], ...
                    %     'parcellation',parcellation, ...
                    %     'surface_name', fsa, ...
                    %     'color_range',[lim_min lim_max], ...
                    %     'label_text',['InterpolationResults_',var,variable,x_var{:}], ...
                    %     'position_colorbar','East')
                    % colorbar_white_centered([lim_min lim_max])

                    % plot_cortical( ...
                    %     parcel_to_surface(Interpolation_Results{idx_to_ENIGMA,['res_',x_var{:}]}, [parcellation,'_',fsa]), ...
                    %     'surface_name', fsa, ...
                    %     'color_range',[lim_min lim_max], ...
                    %     'label_text',['InterpolationResults_',var,variable,x_var{:}])
                    % colorbar_white_centered([lim_min lim_max])

                    % plot_cortical( ...
                    %     parcel_to_surface(Interpolation_Results{idx_to_ENIGMA,['res_',x_var{:}]}, [parcellation,'_',fsa]), ...
                    %     'surface_name', fsa, ...
                    %     'color_range',[-max(abs(Interpolation_Results{:,['res_complete_',x_var{:}]})) max(abs(Interpolation_Results{:,['res_complete_',x_var{:}]}))], ...
                    %     'label_text',['InterpolationResults_',var,variable,x_var{:}])
                    % colorbar_white_centered([lim_min lim_max])

                    % plot_cortical( ...
                    %     parcel_to_surface(Interpolation_Results{idx_to_ENIGMA,['res_complete_',x_var{:}]}, [parcellation,'_',fsa]), ...
                    %     'surface_name', fsa, ...
                    %     'color_range',[lim_min lim_max], ...
                    %     'label_text',['InterpolationResults_',var,variable,x_var{:}], ...
                    %     'position_colorbar','East')
                    % colorbar_white_centered([lim_min lim_max])
                % end
            
            end 
        end
    end
end
