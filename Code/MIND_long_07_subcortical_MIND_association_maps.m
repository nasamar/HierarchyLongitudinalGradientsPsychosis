%% Script to plot brain maps of subcortical MIND associations.

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

parcellation = 'subcortical';

degree_68_FEP = readtable([location,'Code\3. MIND_long\Data\degree\',parcellation,'\COMBATLS_covars\degree_68_FEP.csv'],ReadRowNames=true);
degree_68_CN = readtable([location,'Code\3. MIND_long\Data\degree\',parcellation,'\COMBATLS_covars\degree_68_CN.csv'],ReadRowNames=true);


names_DK = {'L-accumbens', 'L-amygdala', 'L-caudate', 'L-hippocampus', 'L-pallidum', 'L-putamen', 'L-thalamus', 'L-ventricle', ...
    'R-accumbens', 'R-amygdala', 'R-caudate', 'R-hippocampus', 'R-pallidum', 'R-putamen', 'R-thalamus', 'R-ventricle'};
names_DK = replace(replace(names_DK,'R-','rh_'),'L-','lh_');

% Interpolation results
var = 'degrees';
% for period = {'_baseline'}
for period = {'_baseline',''}
    % for cognition = {''}
    for cognition = {'','_BPRS'}

        variable = [period{:},cognition{:}];
        
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

        
        opts = detectImportOptions([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],'ReadVariableNames',true);
        opts = setvartype(opts, strcat('res_',x_vars), 'double'); 
        Interpolation_Results = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],opts);
        

        % poolobj = gcp('nocreate');
        % delete(poolobj);
        % parpool(length(x_vars))
        
        var_partial = mean(degree_68_CN{contains(degree_68_CN.Properties.RowNames,'_001'),:})';
        var_name = 'degree CN';

        lim_max = max(max(Interpolation_Results{:,cellfun(@(x) ['res_complete_' x], x_vars, 'UniformOutput', false)}));
        lim_min = min(min(Interpolation_Results{:,cellfun(@(x) ['res_complete_' x], x_vars, 'UniformOutput', false)}));
        
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
            figure('Position', [488   242   560  200])
               plot_subcortical(Interpolation_Results{:,['res_complete_',x_var{:}]}',...
                    'ventricles','False',...
                    'color_range',[lim_min lim_max], ...
                    'label_text',['InterpolationResults_',var,variable,x_var{:}]);
                colorbar_white_centered([lim_min lim_max])

            
            % figure('Position', [488   242   560  200])
            %    plot_subcortical([Interpolation_Results{:,['res_complete_',x_var{:}]}';Interpolation_Results{:,['res_',x_var{:}]}'],...
            %         'ventricles','False',...
            %         'color_range',[lim_min lim_max], ...
            %         'label_text',['InterpolationResults_',var,variable,x_var{:}]);
            %     colorbar_white_centered([lim_min lim_max])                
                
            % end
            
        end
        
    end

end

