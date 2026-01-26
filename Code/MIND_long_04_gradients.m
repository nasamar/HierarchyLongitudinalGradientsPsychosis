%% Script to compute HC and SSD MIND gradients.

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

clear
close all
location = 'C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\';
type = 'COMBATLS_covars';

residuals = false; % false for raw data

parcellation = 'aparc_500_sym';
% parcellation = 'subcortical';

normalization = 'FEP+CN';
% normalization = 'CN';

edge_metrics = {'edge'};
regional_metrics = {'degree'};
metrics = [edge_metrics,regional_metrics];

for i = 1:length(metrics)
    reading_path = ['Code\3. MIND_long\Data\',metrics{i},'\',parcellation,'\',type,'\'];
    eval([['read_',metrics{i}] ' = reading_path;']);
end


read_edge = ['Code\3. MIND_long\Data\edges\',parcellation,'\',type,'\'];
read_effsizes_edges = ['Code\2. MIND\Data\edges\',parcellation,'\',type,'\'];


for i = 1:length(regional_metrics)
    reading_path_CN = readtable([location,eval(['read_',regional_metrics{i}]),[regional_metrics{i},'_68_CN.csv']],"ReadRowNames",true);
    eval([[regional_metrics{i},'_68_CN'] ' = reading_path_CN;']);
    reading_path_FEP = readtable([location,eval(['read_',regional_metrics{i}]),[regional_metrics{i},'_68_FEP.csv']],"ReadRowNames",true);
    eval([[regional_metrics{i},'_68_FEP'] ' = reading_path_FEP;']);

    % Comment for raw MIND (not residuals)
    if strcmp(regional_metrics{i},'degree') & residuals
        reading_path_CN = readtable([location,eval(['read_',regional_metrics{i}]),[regional_metrics{i},'_68_CN_residuals.csv']],"ReadRowNames",true);
        eval([[regional_metrics{i},'_68_CN'] ' = reading_path_CN;']);
        reading_path_FEP = readtable([location,eval(['read_',regional_metrics{i}]),[regional_metrics{i},'_68_FEP_residuals.csv']],"ReadRowNames",true);
        eval([[regional_metrics{i},'_68_FEP'] ' = reading_path_FEP;']);

    end
end

for i = 1:length(edge_metrics)
    
    load([location,eval(['read_',edge_metrics{i}]),[edge_metrics{i},'_68_CN.mat']])
    load([location,eval(['read_',edge_metrics{i}]),[edge_metrics{i},'_68_FEP.mat']])

    % Comment for raw MIND (not residuals)
    if residuals
        load([location,eval(['read_',edge_metrics{i}]),[edge_metrics{i},'_long_CN_residuals.mat']])
        load([location,eval(['read_',edge_metrics{i}]),[edge_metrics{i},'_long_FEP_residuals.mat']])
    end
    
end

mean_edge_68_FEP = zeros(width(degree_68_FEP));
mean_edge_68_CN = zeros(width(degree_68_CN));
mean_edge_68_FEP_baseline = mean_edge_68_FEP;
mean_edge_68_CN_baseline = mean_edge_68_CN;
edge_68_FEP_baseline = edge_68_FEP(contains(edge_68_FEP(:,1),'001'),:);
edge_68_CN_baseline = edge_68_CN(contains(edge_68_CN(:,1),'001'),:);

for i = 1:length(edge_68_FEP)
    mean_edge_68_FEP = mean_edge_68_FEP + edge_68_FEP{i,2}{:,:};
end
for i = 1:length(edge_68_FEP_baseline)
    mean_edge_68_FEP_baseline = mean_edge_68_FEP_baseline + edge_68_FEP_baseline{i,2}{:,:};
end
for i = 1:length(edge_68_CN)
    mean_edge_68_CN = mean_edge_68_CN + edge_68_CN{i,2}{:,:};
end
for i = 1:length(edge_68_CN_baseline)
    mean_edge_68_CN_baseline = mean_edge_68_CN_baseline + edge_68_CN_baseline{i,2}{:,:};
end
mean_edge_68_FEP = array2table(mean_edge_68_FEP/length(edge_68_FEP));
mean_edge_68_CN = array2table(mean_edge_68_CN/length(edge_68_CN));
mean_edge_68_baseline = (mean_edge_68_FEP_baseline+mean_edge_68_CN_baseline)/(length(edge_68_FEP_baseline)+length(edge_68_CN_baseline));
mean_edge_68_baseline = array2table(mean_edge_68_baseline);
mean_edge_68_FEP_baseline = array2table(mean_edge_68_FEP_baseline/length(edge_68_FEP_baseline));
mean_edge_68_CN_baseline = array2table(mean_edge_68_CN_baseline/length(edge_68_CN_baseline));


variables_embeding = {'mean_edge_68_FEP','mean_edge_68_CN','mean_edge_68_FEP_baseline','mean_edge_68_CN_baseline'};
variables_embeding = {'mean_edge_68_FEP_baseline','mean_edge_68_CN_baseline'};


volumes = readtable([location,'Code\2. MIND\Data\volumes\CorticalMeasuresENIGMA_GrayAvg.csv'],'ReadRowNames',true);
[volumes_names,idx] = sort(volumes(:,1:68).Properties.VariableNames);


if strcmp(parcellation,'aparc_500_sym')
    aparc_TO_500_sym_names = regexp(degree_68_FEP.Properties.VariableNames, '(?<=^[^_]*_)[^_]*(?=_)', 'match', 'once');
    [~,idx] = ismember(aparc_TO_500_sym_names,unique(aparc_TO_500_sym_names));
    
    names_DK = readtable([location,'Molecular\parcellations\aparc_500_sym\500.sym_names.txt'],"ReadVariableNames",false);
    names_DK = cellfun(@(x, y, z) strcat(x,'_',y,'_', z), names_DK{:,1},names_DK{:,2},names_DK{:,3}, 'UniformOutput', false);    

elseif strcmp(parcellation,'subcortical')
    names_DK = degree_68_FEP.Properties.VariableNames;
end

[~,idx_ordered] = ismember(names_DK,degree_68_FEP.Properties.VariableNames);
[~,idx] = ismember(degree_68_FEP.Properties.VariableNames,names_DK);



%% .ANNOT (region) to surface 
if strcmp(parcellation,'aparc_500_sym')
    annot_file = '500.sym.aparc';
    fsa = 'fsa';

    g_lh = gifti(['C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\MATLAB\Libraries\ENIGMA-master\ENIGMA-master\matlab\shared\surfaces\',fsa,'_lh.gii']);
    lh_S.coord = g_lh.vertices';
    lh_S.tri = g_lh.faces;
    
    g_rh = gifti(['C:\Users\usuario\OneDrive - UNIVERSIDAD DE SEVILLA\Natalia\Code\MATLAB\Libraries\ENIGMA-master\ENIGMA-master\matlab\shared\surfaces\',fsa,'_rh.gii']);
    rh_S.coord = g_rh.vertices';
    rh_S.tri = g_rh.faces;
    
    lh_rh_S.coord = [lh_S.coord,rh_S.coord]; 
    lh_rh_S.tri = [lh_S.tri;rh_S.tri]; 


    [lh_vertices, lh_labels, lh_colortable] = read_annotation([location,'Molecular\parcellations\',parcellation,'\lh.',annot_file,'.annot']);
    [rh_vertices, rh_labels, rh_colortable] = read_annotation([location,'Molecular\parcellations\',parcellation,'\rh.',annot_file,'.annot']);
    
    nVert = length(lh_vertices);  
    
    % Inicializar celdas para nombres de región
    lh_region = cell(nVert,1);
    rh_region = cell(nVert,1);
    
    for i = 1:nVert            
        % Mapear etiquetas de vértices a nombres de región (LH)
        idx_vertex = find(lh_colortable.table(:,5) == lh_labels(i));
        if ~isempty(idx_vertex)
            lh_region{i} = lh_colortable.struct_names{idx_vertex};
        else
            lh_region{i} = 'unknown';
        end
    
        % Mapear etiquetas de vértices a nombres de región (RH)
        idx_vertex = find(rh_colortable.table(:,5) == rh_labels(i));
        if ~isempty(idx_vertex)
            rh_region{i} = rh_colortable.struct_names{idx_vertex};
        else
            rh_region{i} = 'unknown';
        end
        % idx_vertex
    end
    
    
    % Crear tablas
    T_lh = table(repmat("lh", nVert, 1), (0:nVert-1)', lh_region, ...
        'VariableNames', {'Hemisphere', 'VertexIndex', 'Region'});
    
    T_rh = table(repmat("rh", nVert, 1), (0:nVert-1)', rh_region, ...
        'VariableNames', {'Hemisphere', 'VertexIndex', 'Region'});
    
    % Concatenar
    T = [T_lh; T_rh];
    
    T = strcat(T{:,1},'_',T{:,3});
    [~,idx_vertex] = ismember(T,names_DK);
    
    
    if strcmp(parcellation,'aparc') & strcmp(fsa,'fsa5') 
        idx_vertex = table2array(readtable([location,'Code\MATLAB\Libraries\ENIGMA-master\ENIGMA-master\matlab\shared\parcellations\aparc_fsa5.csv']));
        idx_vertex(idx_vertex==4) = 0;
        idx_vertex(idx_vertex==39) = 0;
        idx_vertex(idx_vertex > 4) = idx_vertex(idx_vertex > 4) - 1;    
        idx_vertex(idx_vertex > 38) = idx_vertex(idx_vertex > 38) - 1;
    
        writetable(array2table(idx_vertex), [location,'Molecular\parcellations\',parcellation,'\',parcellation,'_fsa5.csv'],'WriteVariableNames',false);
    
    else
        writetable(array2table(idx_vertex), [location,'Molecular\parcellations\',parcellation,'\',parcellation,'_fsa.csv'],'WriteVariableNames',false);
        writetable(array2table(idx_vertex), [location,'Code\MATLAB\Libraries\ENIGMA-master\ENIGMA-master\matlab\shared\parcellations\',parcellation,'_fsa.csv'],'WriteVariableNames',false);
    end
else
    [aseg_vol, M, mr_parms, volsz] = load_mgh([location,'Molecular\parcellations\',parcellation,'aseg']);
    idx_vertex = table2array(readtable([location,'Code\MATLAB\Libraries\ENIGMA-master\ENIGMA-master\matlab\shared\parcellations\aparc_aseg_fsa5_with_sctx.csv']));

end
%%            

% % Parcellation
% h = plot_hemispheres(idx_vertex, {lh_S,rh_S}); % edit process_views to change view 
% colormap(h.handles.figure,lines(68))


% GROUP LEVEL
% if ~strcmp(parcellation,'subcortical')
if false
    for variable_embeding = variables_embeding
        W = table2array(eval(variable_embeding{:}));
        W(isnan(W)) = 0;
    
       % Gradient template
        gm_ref = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');
        if ~contains(normalization,'FEP')
                W_ref = eval('mean_edge_68_CN_baseline');
        elseif ~contains(normalization,'+')
            W_ref = eval(['mean_edge_68_',case_control{end},'_baseline']);
    
        else
            W_ref = eval('mean_edge_68_baseline');
        end
        gm_ref = gm_ref.fit(W_ref{:,:});
    
    
        % Calcular gradientes
        gm = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');    
        gm = gm.fit(W);
    
        [d, Z, transform] = procrustes(gm_ref.gradients{1}(:,1:2), gm.gradients{1}(:,1:2));
    
        gradients_1_2 = gm.gradients{:,:}(idx_ordered,1:2);   

        scree_plot(gm.lambda{1});
        title(variable_embeding,'Interpreter','none')
    
        plot_hemispheres(gradients_1_2,{lh_S,rh_S}, 'parcellation', idx_vertex, 'labeltext',{'Gradient 1','Gradient 2'});  % edit process_views to change view
        % plot_hemispheres([gradient_1',gradient_2'],{lh_S,rh_S}, 'labeltext',{'Gradient 1','Gradient 2'});  % edit process_views to change view
        enigma_colormap('RdBu_r')

        gradient_1 = parcel_to_surface(gradients_1_2(:,1), [parcellation,'_',fsa]);
        gradient_2 = parcel_to_surface(gradients_1_2(:,2),  [parcellation,'_',fsa]);
    
        % figure('Position', [488   242   560  200])
        % plot_cortical(gradient_1, 'surface_name', fsa, 'color_range',[-max(abs(gradient_1)) max(gradient_1)], 'cmap', 'RdBu_r')
        % sgtitle([variable_embeding,' PC1'])
        % 
        % figure('Position', [488   242   560  200])
        % plot_cortical(gradient_2, 'surface_name', fsa, 'color_range',[-max(abs(gradient_2)) max(abs(gradient_2))], 'cmap', 'RdBu_r')
        % sgtitle([variable_embeding,' PC2'])
    
         
        % Half brain
        [h,C] = gradient_in_euclidean(gradients_1_2,{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
        % [h,C] = gradient_in_euclidean(gradients,{lh_S,rh_S},[idx_vertex(1:length(idx_vertex)/2),idx_vertex(length(idx_vertex)/2+1:end)]);
        sgtitle(variable_embeding,'Interpreter','none')
        % xlim([-0.2 0.2])
        % ylim([-0.2 0.2])
    
        [h,C] = gradient_in_euclidean(gm_ref.gradients{1}(idx_ordered,1:2),{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
        sgtitle(normalization,'Interpreter','none')
    
        [h,C] = gradient_in_euclidean(Z(idx_ordered,:),{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
        sgtitle([variable_embeding,' normalized'],'Interpreter','none')

        % Full brain
        % [h,C] = gradient_in_euclidean(gradients,lh_rh_S,idx_vertex);
        % sgtitle(variable_embeding,'Interpreter','none')

        figure;
        scatter(Z(idx_ordered,1),Z(idx_ordered,2),'filled')
        xlabel('G1')
        ylabel('G2')
        hold on
        legend({'FEP baseline','CN baseline'})
     

        
        % close all
    end
end

if residuals
    residual_or_raw = 'residuals';
else
    residual_or_raw = 'raw';
end
    
% INDIVIDUAL LEVEL
variables_embeding = {'edge_68_FEP','edge_68_CN'};
% variables_embeding = {'edge_68_CN'};

% if ~strcmp(parcellation,'subcortical')
if false
% if ~exist([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G1_CN.csv'])
    for variable_embeding = variables_embeding
        variable = eval(variable_embeding{:});
        case_control = split(variable_embeding{:},'_'); 

        gradient_G1 = [];
        gradient_G2 = [];
        SD = [];
        num_inside = [];

        gradient_G1_normalized = [];
        gradient_G2_normalized = [];
        SD_normalized = [];
        num_inside_normalized = [];
        range_G1_normalized = [];
        range_G2_normalized = [];

        std_G1_normalized = [];
        std_G2_normalized = [];

        gm_ref = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');
        if ~contains(normalization,'FEP')
            W_ref = eval('mean_edge_68_CN_baseline');
        elseif ~contains(normalization,'+')
            W_ref = eval(['mean_edge_68_',case_control{end},'_baseline']);

        else
            W_ref = eval('mean_edge_68_baseline');
        end
        gm_ref = gm_ref.fit(W_ref{:,:});


        for i = 1:length(variable)
            
            W = table2array(variable{i,2});

            % Calcular gradientes
            gm = GradientMaps('approach', 'dm', 'kernel', 'normalized_angle');
            gm = gm.fit(W);
            
            % writetable(array2table(gm.gradients{1}(:,1:2)',"RowNames",{'G1','G2'},"VariableNames",degree_68_FEP.Properties.VariableNames),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients\',parcellation,'\',type,'\',residual_or_raw,'\',variable{i,1},'.csv'],"WriteRowNames",true)

            [d, Z, transform] = procrustes(gm_ref.gradients{1}(:,1:2), gm.gradients{1}(:,1:2));
            writetable(array2table(Z',"RowNames",{'G1','G2'},"VariableNames",degree_68_FEP.Properties.VariableNames),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\',variable{i,1},'.csv'],"WriteRowNames",true)

            gradient_G1(i,:) = gm.gradients{1}(:,1)';
            eval([['gradient_G1_',case_control{end}] ' = gradient_G1;']);

            gradient_G2(i,:) = gm.gradients{1}(:,2)';
            eval([['gradient_G2_',case_control{end}] ' = gradient_G2;']);

            SD(i,:) = sqrt(mean((gm.gradients{1}(:,1) - mean(gm.gradients{1}(:,1))).^2 + (gm.gradients{1}(:,2) - mean(gm.gradients{1}(:,2))).^2)); % standard distance
            eval([['SD_',case_control{end}] ' = SD;']);

            num_inside(i,:) = sum((gm.gradients{1}(:,1) - mean(gm.gradients{1}(:,1))).^2 + (gm.gradients{1}(:,2) - mean(gm.gradients{1}(:,2))).^2 <= SD(i,:)^2);
            eval([['num_inside_',case_control{end}] ' = num_inside;']);


            gradient_G1_normalized(i,:) = Z(:,1)';

            gradient_G2_normalized(i,:) = Z(:,2)';

            SD_normalized(i,:) = sqrt(mean((Z(:,1) - mean(Z(:,1))).^2 + (Z(:,2) - mean(Z(:,2))).^2)); % standard distance
            eval([['SD_normalized_',case_control{end}] ' = SD_normalized;']);

            num_inside_normalized(i,:) = sum((Z(:,1) - mean(Z(:,1))).^2 + (Z(:,2) - mean(Z(:,2))).^2 <= SD(i,:)^2);
            eval([['num_inside_normalized_',case_control{end}] ' = num_inside_normalized;']);

            range_G1_normalized(i,:) = range(Z(:,1));
            eval([['range_G1_normalized_',case_control{end}] ' = range_G1_normalized;']);

            range_G2_normalized(i,:) = range(Z(:,2));
            eval([['range_G2_normalized_',case_control{end}] ' = range_G2_normalized;']);

            std_G1_normalized(i,:) = std(Z(:,1));
            eval([['std_G1_normalized_',case_control{end}] ' = std_G1_normalized;']);

            std_G2_normalized(i,:) = std(Z(:,2));
            eval([['std_G2_normalized_',case_control{end}] ' = std_G2_normalized;']);
            
            if i == length(variable)
                gradient_G1_normalized = gradient_G1_normalized(any(imag(gradient_G1_normalized) == 0, 2),:);
                gradient_G2_normalized = gradient_G2_normalized(any(imag(gradient_G2_normalized) == 0, 2),:);
                variable = variable(any(imag(gradient_G2_normalized) == 0, 2),:);

            end
            eval([['gradient_G1_normalized_',case_control{end}] ' = gradient_G1_normalized;']);
            eval([['gradient_G2_normalized_',case_control{end}] ' = gradient_G2_normalized;']);


            % Gradients
            gradientes = gm.gradients{1}(:,1:2);
            

            % % PLOTS
            gradients_1_2 = gm.gradients{:,:}(idx_ordered,1:2);


            % scree_plot(gm.lambda{1});
            % title([variable_embeding,' ',variable{i,1}],'Interpreter','none')


            % plot_hemispheres(gradients_1_2,{lh_S,rh_S}, 'parcellation', idx_vertex, 'labeltext',{'Gradient 1','Gradient 2'});  % edit process_views to change view
            % enigma_colormap('RdBu_r')

            % gradient_1 = parcel_to_surface(gradients_1_2(:,1),  [parcellation,'_',fsa]);
            % gradient_2 = parcel_to_surface(gradients_1_2(:,2), [parcellation,'_',fsa]);


            % figure;
            % plot_cortical(gradient_1, 'surface_name', fsa, 'color_range',[min(gradient_1) max(gradient_1)], 'cmap', 'RdBu_r')
            % sgtitle([variable_embeding,' ',variable{i,1},' PC1'])
            % 
            % figure;
            % plot_cortical(gradient_2, 'surface_name', fsa, 'color_range',[min(gradient_2) max(gradient_2)], 'cmap', 'RdBu_r')
            % sgtitle([variable_embeding,' ',variable{i,1},' PC2'])


            % % Half brain
            % [h,C] = gradient_in_euclidean(gradients_1_2,{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
            % % [h,C] = gradient_in_euclidean(gradients_1_2,{lh_S,rh_S},[idx_vertex(1:length(idx_vertex)/2),idx_vertex(length(idx_vertex)/2+1:end)]);
            % sgtitle([variable_embeding,' ',variable{i,1}],'Interpreter','none')
            % % xlim([-0.2 0.2])
            % % ylim([-0.2 0.2])
            % 
            % [h,C] = gradient_in_euclidean(gm_ref.gradients{1}(idx_ordered,1:2),{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
            % sgtitle([variable_embeding,' ',variable{i,1}],'Interpreter','none')
            % 
            % [h,C] = gradient_in_euclidean(Z(idx_ordered,:),{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
            % sgtitle([variable_embeding,' ',variable{i,1}],'Interpreter','none')
            % 
            % Full brain
            % [h,C] = gradient_in_euclidean(gradients_1_2,lh_rh_S,idx_vertex);
            % sgtitle([variable_embeding,' ',variable{i,1}],'Interpreter','none')

            
            % [h,C] = gradient_in_euclidean(Z(idx_ordered,:)); % edit gradient_in_euclidean to change view
            % sgtitle([variable_embeding,' ',variable{i,1}],'Interpreter','none') 
            % saveas(gcf,[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\plots\',case_control{end},'\',variable{i,1},'.png'])
            
            

            % Plot circle
            % theta = linspace(0, 2*pi, 100);
            % xcirc = mean(Z(:,1)) + SD_normalized * cos(theta);
            % ycirc = mean(Z(:,2)) + SD_normalized * sin(theta);
            % hold on
            % plot(xcirc, ycirc, 'r-', 'LineWidth', 2)
            % plot(mean(Z(:,1)), mean(Z(:,2)), 'rx', 'MarkerSize', 10, 'LineWidth', 2)  % centroid
            % 
            % 
            % close all

        end

        writetable(array2table(eval(['gradient_G1_normalized_',case_control{end}]),"RowNames",variable(:,1),"VariableNames",degree_68_FEP.Properties.VariableNames),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G1_',case_control{end},'.csv'],"WriteRowNames",true)
        writetable(array2table(eval(['gradient_G2_normalized_',case_control{end}]),"RowNames",variable(:,1),"VariableNames",degree_68_FEP.Properties.VariableNames),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G2_',case_control{end},'.csv'],"WriteRowNames",true)
        % writetable(array2table(eval(['SD_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\SD_',case_control{end},'.csv'],"WriteRowNames",true)
        % writetable(array2table(eval(['num_inside_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\num_inside_',case_control{end},'.csv'],"WriteRowNames",true)
        % 
        % writetable(array2table(eval(['range_G1_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G1_',case_control{end},'.csv'],"WriteRowNames",true)
        % writetable(array2table(eval(['range_G2_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G2_',case_control{end},'.csv'],"WriteRowNames",true)
        % 
        % writetable(array2table(eval(['std_G1_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G1_',case_control{end},'.csv'],"WriteRowNames",true)
        % writetable(array2table(eval(['std_G2_normalized_',case_control{end}]),"RowNames",variable(:,1)),[location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G2_',case_control{end},'.csv'],"WriteRowNames",true)

    end
     
else
    gradients_G1_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G1_FEP.csv'],"ReadRowNames",true);
    gradients_G1_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G1_CN.csv'],"ReadRowNames",true);
    gradients_G2_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G2_FEP.csv'],"ReadRowNames",true);
    gradients_G2_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\gradients_G2_CN.csv'],"ReadRowNames",true);

    SD_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\SD_FEP.csv'],"ReadRowNames",true);
    SD_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\SD_CN.csv'],"ReadRowNames",true);

    % range_G1_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G1_FEP.csv'],"ReadRowNames",true);
    % range_G1_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G1_CN.csv'],"ReadRowNames",true);
    % range_G2_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G2_FEP.csv'],"ReadRowNames",true);
    % range_G2_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\range_G2_CN.csv'],"ReadRowNames",true);
    % 
    % std_G1_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G1_FEP.csv'],"ReadRowNames",true);
    % std_G1_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G1_CN.csv'],"ReadRowNames",true);
    % std_G2_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G2_FEP.csv'],"ReadRowNames",true);
    % std_G2_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\std_G2_CN.csv'],"ReadRowNames",true);

    
    % num_inside_FEP = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\num_inside_FEP.csv'],"ReadRowNames",true);
    % num_inside_CN = readtable([location,'Code\MATLAB\Connectivity\Longitudinal\gradients_normalized\',parcellation,'\',type,'\',normalization,'\',residual_or_raw,'\num_inside_CN.csv'],"ReadRowNames",true);


end

opts = detectImportOptions([location,'\datasets\PAFIP\Covariates_complete.csv'],ReadRowNames=true);
opts = setvartype(opts, 'Machine_Teslas', 'string'); 

covariates = readtable([location,'\datasets\PAFIP\Covariates_complete.csv'],opts);

covariates_FEP = covariates(edge_68_FEP(:,1),:);
covariates_CN = covariates(edge_68_CN(:,1),:);



if strcmp(parcellation,'aparc_500_sym')

    % plot_cortical(parcel_to_surface(mean(degree_68_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND FEP')
    plot_cortical(parcel_to_surface(mean(degree_68_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[0.1 max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND FEP')
    colorbar_white_centered([min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))])

    % plot_cortical(parcel_to_surface(mean(degree_68_CN{covariates_CN.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND CN')
    plot_cortical(parcel_to_surface(mean(degree_68_CN{covariates_CN.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[0.1 max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND CN')
    colorbar_white_centered([min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))])
    
    % plot_cortical(parcel_to_surface(mean(degree_68_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa],median(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))), 'surface_name', fsa, 'color_range',[min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND FEP')
    % colorbar_white_centered([min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))],median(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))))
    % 
    % plot_cortical(parcel_to_surface(mean(degree_68_CN{covariates_CN.Assessment==1,idx_ordered}),[parcellation,'_',fsa],median(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))), 'surface_name', fsa, 'color_range',[min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','MIND CN')
    % colorbar_white_centered([min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))],median(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))))
    

    % plot_hemispheres(parcel_to_surface(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa])',{lh_S,rh_S},'labeltext',{'Gradient 1'});  % edit process_views to change view
    % enigma_colormap('RdBu_r')

    plot_cortical(parcel_to_surface(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[-max(abs(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','G1 FEP')
    colorbar_white_centered([-max(abs(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,:})))])

    plot_cortical(parcel_to_surface(mean(gradients_G2_FEP{covariates_FEP.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[-max(abs(mean(gradients_G2_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(gradients_G2_FEP{covariates_FEP.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','G2 FEP')
    colorbar_white_centered([-max(abs(mean(gradients_G2_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(gradients_G2_FEP{covariates_FEP.Assessment==1,:})))])

    plot_cortical(parcel_to_surface(mean(gradients_G1_CN{covariates_CN.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[-max(abs(mean(gradients_G1_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(gradients_G1_CN{covariates_CN.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','G1 CN')
    colorbar_white_centered([-max(abs(mean(gradients_G1_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(gradients_G1_CN{covariates_CN.Assessment==1,:})))])

    plot_cortical(parcel_to_surface(mean(gradients_G2_CN{covariates_CN.Assessment==1,idx_ordered}),[parcellation,'_',fsa]), 'surface_name', fsa, 'color_range',[-max(abs(mean(gradients_G2_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(gradients_G2_CN{covariates_CN.Assessment==1,:})))], 'cmap', 'RdBu_r','label_text','G2 CN')
    colorbar_white_centered([-max(abs(mean(gradients_G2_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(gradients_G2_CN{covariates_CN.Assessment==1,:})))])


    [h,C] = gradient_in_euclidean([mean(gradients_G1_FEP{covariates_FEP.Assessment==1,idx_ordered})',mean(gradients_G2_FEP{covariates_FEP.Assessment==1,idx_ordered})'],{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
    sgtitle(['FEP baseline normalized'],'Interpreter','none')

    [h,C] = gradient_in_euclidean([mean(gradients_G1_CN{covariates_CN.Assessment==1,idx_ordered})',mean(gradients_G2_CN{covariates_CN.Assessment==1,idx_ordered})'],{lh_S,rh_S},idx_vertex); % edit gradient_in_euclidean to change view
    sgtitle(['CN baseline normalized'],'Interpreter','none')
    

    figure;
    scatter(mean(gradients_G1_FEP{covariates_FEP.Assessment==1,idx_ordered})',mean(gradients_G2_FEP{covariates_FEP.Assessment==1,idx_ordered})','filled')
    hold on
    scatter(mean(gradients_G1_CN{covariates_CN.Assessment==1,idx_ordered})',mean(gradients_G2_CN{covariates_CN.Assessment==1,idx_ordered})','filled')
    xlabel('G1')
    ylabel('G2')
    legend({'FEP baseline','CN baseline'})
    
    figure;
    imagesc(gradients_G1_CN{:,:})
    title('G1 CN')
    colorbar
    % 
    % figure;
    % imagesc(gradients_G1_FEP{:,:})
    % title('G1 FEP')
    % colorbar
    % 
    % figure;
    % imagesc(gradients_G2_CN{:,:})
    % title('G2 CN')
    % colorbar
    % 
    % figure;
    % imagesc(gradients_G2_FEP{:,:})
    % title('G2 FEP')
    % colorbar
    
    % figure;
    % imagesc(gradients_G1_CN{covariates_CN.Assessment==1,:})
    % title('G1 CN baseline')
    % colorbar
    % 
    % figure;
    % imagesc(gradients_G1_FEP{covariates_FEP.Assessment==1,:})
    % title('G1 FEP baseline')
    % colorbar
    % 
    % figure;
    % imagesc(gradients_G2_CN{covariates_CN.Assessment==1,:})
    % title('G2 CN baseline')
    % colorbar
    % 
    % figure;
    % imagesc(gradients_G2_FEP{covariates_FEP.Assessment==1,:})
    % title('G2 FEP baseline')
    % colorbar
    
    % figure;
    % scatter(mean(gradients_G1_CN{covariates_CN.Assessment==1,:})', mean(gradients_G1_FEP{covariates_FEP.Assessment==1,:})','filled');
    % lsline
    % xlabel('G1 CN baseline')
    % ylabel('G1 FEP baseline')
    % 
    % figure;
    % scatter(mean(gradients_G2_CN{covariates_CN.Assessment==1,:})', mean(gradients_G2_FEP{covariates_FEP.Assessment==1,:})','filled');
    % lsline
    % xlabel('G2 CN baseline')
    % ylabel('G2 FEP baseline')
    
    clims = [min(min([corr(gradients_G1_CN{:,:},gradients_G2_CN{:,:}),corr(gradients_G1_FEP{:,:},gradients_G2_FEP{:,:})])), 
        max(max([corr(gradients_G1_CN{:,:},gradients_G2_CN{:,:}),corr(gradients_G1_FEP{:,:},gradients_G2_FEP{:,:})]))];
    
    figure;
    imagesc(corr(gradients_G1_CN{:,:},gradients_G2_CN{:,:}),clims)
    title('CN')
    colorbar
    
    % figure;
    % imagesc(corr(gradients_G1_FEP{:,:},gradients_G2_FEP{:,:}),clims)
    % title('FEP')
    % colorbar
    % 
    % clims = [min([min(min(corr(gradients_G1_CN{:,:}',gradients_G2_CN{:,:}'))),min(min(corr(gradients_G1_FEP{:,:}',gradients_G2_FEP{:,:}')))]), 
    %     max([max(max(corr(gradients_G1_CN{:,:}',gradients_G2_CN{:,:}'))),max(max(corr(gradients_G1_FEP{:,:}',gradients_G2_FEP{:,:}')))])];
    % 
    % figure;
    % imagesc(corr(gradients_G1_CN{:,:}',gradients_G2_CN{:,:}'),[-0.88 0.51])
    % title('CN')
    % colorbar
    % 
    % figure;
    % imagesc(corr(gradients_G1_FEP{:,:}',gradients_G2_FEP{:,:}'),[-0.88 0.51])
    % title('FEP')
    % colorbar
    
    
    % clims = [min(min([corr(gradients_G1_CN{covariates_CN.Assessment==1,:},gradients_G2_CN{covariates_CN.Assessment==1,:}),corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:},gradients_G2_FEP{covariates_FEP.Assessment==1,:})])), 
    %     max(max([corr(gradients_G1_CN{covariates_CN.Assessment==1,:},gradients_G2_CN{covariates_CN.Assessment==1,:}),corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:},gradients_G2_FEP{covariates_FEP.Assessment==1,:})]))];
    % 
    % figure;
    % imagesc(corr(gradients_G1_CN{covariates_CN.Assessment==1,:},gradients_G2_CN{covariates_CN.Assessment==1,:}),[-0.64 0.66])
    % title('CN baseline')
    % colorbar
    % 
    % figure;
    % imagesc(corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:},gradients_G2_FEP{covariates_FEP.Assessment==1,:}),[-0.64 0.66])
    % title('FEP baseline')
    % colorbar
    % 
    % clims = [min([min(min(corr(gradients_G1_CN{covariates_CN.Assessment==1,:}',gradients_G2_CN{covariates_CN.Assessment==1,:}'))),min(min(corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:}',gradients_G2_FEP{covariates_FEP.Assessment==1,:}')))]), 
    %     max([max(max(corr(gradients_G1_CN{covariates_CN.Assessment==1,:}',gradients_G2_CN{covariates_CN.Assessment==1,:}'))),max(max(corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:}',gradients_G2_FEP{covariates_FEP.Assessment==1,:}')))])];
    % 
    % figure;
    % imagesc(corr(gradients_G1_CN{covariates_CN.Assessment==1,:}',gradients_G2_CN{covariates_CN.Assessment==1,:}'),[-0.75 0.39])
    % title('CN baseline')
    % colorbar
    % 
    % figure;
    % imagesc(corr(gradients_G1_FEP{covariates_FEP.Assessment==1,:}',gradients_G2_FEP{covariates_FEP.Assessment==1,:}'),[-0.75 0.39])
    % title('FEP baseline')
    % colorbar


elseif strcmp(parcellation,'subcortical')
    figure('Position', [488   242   560  200])
    plot_subcortical(mean(degree_68_FEP{covariates_FEP.Assessment==1,idx_ordered}), ...
        'ventricles','False',...
        'color_range',[min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))], ...
        'cmap', 'RdBu_r', ...
        'label_text','MIND FEP')                       
    colorbar_white_centered([min(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:}))) max(abs(mean(degree_68_FEP{covariates_FEP.Assessment==1,:})))])

    figure('Position', [488   242   560  200])
    plot_subcortical(mean(degree_68_CN{covariates_CN.Assessment==1,idx_ordered}), ...
        'ventricles','False',...
        'color_range',[min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))], ...
        'cmap', 'RdBu_r', ...
        'label_text','MIND CN')                       
    colorbar_white_centered([min(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:}))) max(abs(mean(degree_68_CN{covariates_CN.Assessment==1,:})))])

end


