%%
clearvars
clc

%% Models AGORA

Bacteroides_thetaiotaomicron_VPI_5482_AGORA  = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482_AGORA');
Bifidobacterium_animalis_lactis_PT33_AGORA = readCbModel('Bifidobacterium_animalis_lactis_PT33_AGORA.mat');
Clostridium_innocuum_HGF2_AGORA = readCbModel('Clostridium_innocuum_HGF2_AGORA.mat');
Clostridium_sp_7_2_43FAA_AGORA = readCbModel('Clostridium_sp_7_2_43FAA_AGORA.mat');
Clostridium_sp_M62_1_AGORA = readCbModel('Clostridium_sp_M62_1_AGORA.mat');
Clostridium_symbiosum_WAL_14673_AGORA = readCbModel('Clostridium_symbiosum_WAL_14673_AGORA.mat');
Lacticaseibacillus_paracasei_M38_AGORA = readCbModel('Lacticaseibacillus_paracasei_M38_AGORA.mat');

%% Models gapseq

Bacteroides_thetaiotaomicron_VPI_5482_gapseq = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482_gapseq');
Bifidobacterium_animalis_subsp_lactis_PT33_gapseq = readCbModel('Bifidobacterium_animalis_subsp_lactis_PT33_gapseq.mat');
Clostridium_innocuum_HFG2_gapseq = readCbModel('Clostridium_innocuum_HFG2_gapseq.mat');
Clostridium_sp_7_2_43FAA_gapseq = readCbModel('Clostridium_sp_7_2_43FAA_gapseq.mat');
Clostridium_sp_M62_gapseq = readCbModel('Clostridium_sp_M62_gapseq.mat');
Clostridium_symbiosum_WAL_14673_gapseq = readCbModel('Clostridium_symbiosum_WAL_14673_gapseq.mat');
Lacticaseibacillus_paracasei_M38_gapseq = readCbModel('Lacticaseibacillus_paracasei_M38_gapseq.mat');

%%
% List of AGORA models
models_AGORA = {Bacteroides_thetaiotaomicron_VPI_5482_AGORA, Bifidobacterium_animalis_lactis_PT33_AGORA, Clostridium_innocuum_HGF2_AGORA, ...
    Clostridium_sp_7_2_43FAA_AGORA, Clostridium_sp_M62_1_AGORA, Clostridium_symbiosum_WAL_14673_AGORA, Lacticaseibacillus_paracasei_M38_AGORA};

% List of gapseq models
models_gapseq = {Bacteroides_thetaiotaomicron_VPI_5482_gapseq, Bifidobacterium_animalis_subsp_lactis_PT33_gapseq, Clostridium_innocuum_HFG2_gapseq, ...
    Clostridium_sp_7_2_43FAA_gapseq, Clostridium_sp_M62_gapseq, Clostridium_symbiosum_WAL_14673_gapseq, Lacticaseibacillus_paracasei_M38_gapseq};

% Data extraction
num_reactions_AGORA = cellfun(@(x) length(x.rxns), models_AGORA);
num_metabolites_AGORA = cellfun(@(x) length(x.mets), models_AGORA);
num_genes_AGORA = cellfun(@(x) length(x.genes), models_AGORA);

num_reactions_gapseq = cellfun(@(x) length(x.rxns), models_gapseq);
num_metabolites_gapseq = cellfun(@(x) length(x.mets), models_gapseq);
num_genes_gapseq = cellfun(@(x) length(x.genes), models_gapseq);

%% Graph metabolites, reactions and genes both sets of models

% Bacteria names for the x-axis
bacteria_names = {'B. thetaiotaomicron', 'B. animalis lactis', 'C. innocuum', 'C. sp 7_2_43FAA', 'C. sp M62_1', 'C. symbiosum', 'L. paracasei'};

% Comparison graphs
figure;

% Number of reactions
subplot(3,1,1);
bar([num_reactions_AGORA; num_reactions_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45);
legend('AGORA', 'gapseq');
ylabel('Número de reacciones');
title('Comparación entre modelos AGORA y gapseq');

% Number of metabolites
subplot(3,1,2);
bar([num_metabolites_AGORA; num_metabolites_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45);
ylabel('Número de metabolitos');

% Number of genes
subplot(3,1,3);
bar([num_genes_AGORA; num_genes_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45);
ylabel('Número de genes');

% Adjust the size of the figure to better visualize the labels
set(gcf, 'Position', [100, 100, 800, 600]);

%% Metabolite formula in common AGORA and gapseq

% List of AGORA models
models_AGORA_list = {Bacteroides_thetaiotaomicron_VPI_5482_AGORA, Bifidobacterium_animalis_lactis_PT33_AGORA, ...
    Clostridium_innocuum_HGF2_AGORA, Clostridium_sp_7_2_43FAA_AGORA, Clostridium_sp_M62_1_AGORA, ...
    Clostridium_symbiosum_WAL_14673_AGORA, Lacticaseibacillus_paracasei_M38_AGORA};

% List of gapseq models
models_gapseq_list = {Bacteroides_thetaiotaomicron_VPI_5482_gapseq, Bifidobacterium_animalis_subsp_lactis_PT33_gapseq, ...
    Clostridium_innocuum_HFG2_gapseq, Clostridium_sp_7_2_43FAA_gapseq, Clostridium_sp_M62_gapseq, ...
    Clostridium_symbiosum_WAL_14673_gapseq, Lacticaseibacillus_paracasei_M38_gapseq};

% Bacteria names for labels
bacteria_names = {'Bacteroides_thetaiotaomicron_VPI_5482', 'Bifidobacterium_animalis_lactis_PT33', ...
    'Clostridium_innocuum_HGF2', 'Clostridium_sp_7_2_43FAA', 'Clostridium_sp_M62_1', ...
    'Clostridium_symbiosum_WAL_14673', 'Lacticaseibacillus_paracasei_M38'};

% Number of models
num_models = length(models_AGORA_list);

% Preallocation of a vector for storing the number of formulas in common
common_formulas_count = zeros(1, num_models);

for i = 1:num_models
    % Extracting metabolic formulas from AGORA models
    formulas_AGORA = models_AGORA_list{i}.metFormulas;

    % Extracting metabolic formulas from gapseq models
    formulas_gapseq = models_gapseq_list{i}.metFormulas;

    % Find formulas in common
    common_formulas = intersect(formulas_AGORA, formulas_gapseq);

    % Store the number of formulas in common
    common_formulas_count(i) = length(common_formulas);
end

% Bar chart showing the number of formulas in common
figure;
bar(common_formulas_count);
set(gca, 'XTick', 1:num_models, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45);
ylabel('Número de fórmulas en común');
title('Fórmulas metabólicas en común entre modelos AGORA y gapseq');

%% 4 subplots for final figure
% Bacteria names for the x-axis
bacteria_names = {'B. thetaiotaomicron VPI 5482', 'B. animalis lactis PT33', 'C. innocuum HGF2', 'C. tertium 7-2-43FAA', 'C. saccharolyticum M62-1', 'C. symbiosum WAL 14673', 'L. paracasei M38'};

% Comparison graphs
figure;

fontSize = 14;  % Define a font size for all axes
y_limit = [0 2500]; % Y-axis limits

% Number of reactions
subplot(2,2,1);
bar([num_reactions_AGORA; num_reactions_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45, 'FontSize', fontSize, 'FontWeight', 'bold');
ylim(y_limit);
legend('AGORA', 'gapseq');
ylabel('Number of reactions', 'FontSize', fontSize, 'FontWeight', 'bold');
title('Reaction comparison', 'FontSize', fontSize);

% Number of metabolites
subplot(2,2,2);
bar([num_metabolites_AGORA; num_metabolites_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45, 'FontSize', fontSize, 'FontWeight', 'bold');
ylim(y_limit);
ylabel('Number of metabolites', 'FontSize', fontSize, 'FontWeight', 'bold');
title('Metabolite comparison', 'FontSize', fontSize);

% Number of genes
subplot(2,2,3);
bar([num_genes_AGORA; num_genes_gapseq]');
set(gca, 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45, 'FontSize', fontSize, 'FontWeight', 'bold');
ylim(y_limit);
ylabel('Number of genes', 'FontSize', fontSize, 'FontWeight', 'bold');
title('Gene comparison', 'FontSize', fontSize);

% Bar chart showing the number of formulas in common
subplot(2,2,4);
bar(common_formulas_count, 'FaceColor', [0 0.6 0]);
set(gca, 'XTick', 1:length(bacteria_names), 'XTickLabel', bacteria_names, 'XTickLabelRotation', 45, 'FontSize', fontSize, 'FontWeight', 'bold');
ylim(y_limit);
ylabel('Number of formulas', 'FontSize', fontSize, 'FontWeight', 'bold');
title('Chemical formulas of common metabolites comparison', 'FontSize', fontSize);

% Adjust the size of the figure to better visualize the labels.
set(gcf, 'Position', [100, 100, 900, 700]);






