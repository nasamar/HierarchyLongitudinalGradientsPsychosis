%% Script to plot brain maps of subcortical MIND associations.

% Copyright (C) 2026 University of Seville

% Written by Natalia García San Martín (ngarcia1@us.es)

% This file is part of Hierarchy Longitudinal MIND Gradients Psychosis toolkit.
%
% Hierarchy Longitudinal MIND Gradients Psychosis toolkit is free software: 
% you can redistribute it and/or modify it under the terms of the 
% GNU General Public License as published by the Free Software Foundation, 
% either version 3 of the License, or (at your option) any later version.
%
% Hierarchy Longitudinal MIND Gradients Psychosis toolkit is distributed in the hope that 
% it will be useful, but WITHOUT ANY WARRANTY; without even the implied 
% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Hierarchy Longitudinal MIND Gradients Psychosis toolkit. If not, see 
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
% for period = {''}
for period = {'_baseline',''}
    % for cognition = {'_BPRS'}
    for cognition = {'','_BPRS'}

        variable = [period{:},cognition{:}];
        
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
            x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler'};
        else
            x_vars = {['',var],'Treatment_Time',[var,'_Treatment_Time'],'CPZ_equivalent','CPZ_equivalent_Treatment_Time'};
            x_vars = {'Age_inclusion','Sex1','eTIV','mean_euler'};            
        end

        
        opts = detectImportOptions([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],'ReadVariableNames',true);
        opts = setvartype(opts, strcat('res_',x_vars), 'double'); 
        Interpolation_Results = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\',parcellation,'\InterpolationResults_',var,variable,'.csv'],opts);
        

        Interpolation_Results_degrees = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\degrees\aparc_500_sym\InterpolationResults_degrees',period{:},cognition{:},'.csv'],opts);
        Interpolation_Results_gradients_G1 = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\aparc_500_sym\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_gradients_G1',period{:},cognition{:},'.csv'],opts);
        Interpolation_Results_gradients_G2 = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\aparc_500_sym\COMBATLS_covars',normalization,'\',residual_or_raw,'\InterpolationResults_gradients_G2',period{:},cognition{:},'.csv'],opts);

        lim_max_degrees = max(max(Interpolation_Results_degrees{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'degrees'), 'UniformOutput', false)}));
        lim_min_degrees = min(min(Interpolation_Results_degrees{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'degrees'), 'UniformOutput', false)}));
        lim_max_G1 = max(max(Interpolation_Results_gradients_G1{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'gradients'), 'UniformOutput', false)}));
        lim_min_G1 = min(min(Interpolation_Results_gradients_G1{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'gradients'), 'UniformOutput', false)}));
        lim_max_G2 = max(max(Interpolation_Results_gradients_G2{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'gradients'), 'UniformOutput', false)}));
        lim_min_G2 = min(min(Interpolation_Results_gradients_G2{:,cellfun(@(x) ['res_complete_' x], replace(x_vars,var,'gradients'), 'UniformOutput', false)}));


        if strcmp(period,'_baseline') | contains(cognition,{'_BPRS'})  
            lim_min = min([lim_min_degrees,lim_min_G1,lim_min_G2]);
            lim_max = max([lim_max_degrees,lim_max_G1,lim_max_G2]);
            
        elseif contains(var,{'degrees'})
            lim_max = lim_max_degrees;
            lim_min = lim_min_degrees;
        else
            lim_max = max([lim_max_G1,lim_max_G2]);
            lim_min = min([lim_min_G1,lim_min_G2]);            
        end

        
        for i = 1:length(x_vars)

            x_var = x_vars(i);

            % Interpolation_Results{isnan(Interpolation_Results{:,['res_',x_var{:}]}),['res_',x_var{:}]} = 0;

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
                
            
        end
        
    end

end
