%% refinement of VMH models for simulation on inulin
% Change working directory if necessary

% initCobraToolbox(false)
changeCobraSolver('gurobi')

%% AGORA v1.3 models downloaded from VMH
Bacteroides_thetaiotaomicron_VPI_5482=readCbModel('Bacteroides_thetaiotaomicron_VPI_5482');
Bifidobacterium_animalis_lactis_AD011=readCbModel('Bifidobacterium_animalis_lactis_AD011');
Clostridium_innocuum_2959=readCbModel('Clostridium_innocuum_2959');
Clostridium_sp_7_2_43FAA=readCbModel('Clostridium_sp_7_2_43FAA');
Clostridium_sp_M62_1=readCbModel('Clostridium_sp_M62_1');
Clostridium_symbiosum_WAL_14673=readCbModel('Clostridium_symbiosum_WAL_14673');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302=readCbModel('Lactobacillus_paracasei_subsp_paracasei_ATCC_25302');

%% Adding fructose consumption reactions to B. animalis subsp. lactis 
% surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "fru[e]");
surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "fru[c]");

surfNet(Bifidobacterium_animalis_lactis_AD011, "fru[e]");

% Perform the same search for the M62 model
surfNet(Clostridium_innocuum_2959, "fru[e]");
surfNet(Clostridium_symbiosum_WAL_14673, "fru[e]");

surfNet(Clostridium_innocuum_2959, "fru[c]");
surfNet(Clostridium_sp_M62_1, "fru[c]");
surfNet(Clostridium_sp_M62_1, "f6p[c]");

surfNet(Clostridium_innocuum_2959, "sucr[e]");
surfNet(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, "sucr[e]");
surfNet(Clostridium_sp_M62_1, "sucr[e]");

surfNet(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, "fru[e]");

%% Addition of reactions according to experimental results (TLC)
% B. animalis subsp. lactis AD011:
Bifidobacterium_animalis_lactis_AD011 = addReaction(Bifidobacterium_animalis_lactis_AD011, 'EX_fru(e)',...
'reactionFormula', 'fru[e]   <=>');
Bifidobacterium_animalis_lactis_AD011 = addReaction(Bifidobacterium_animalis_lactis_AD011, 'FRUpts',...
'fru[e] + pep[c]   ->   f1p[c] + pyr[c]');
Bifidobacterium_animalis_lactis_AD011 = addReaction(Bifidobacterium_animalis_lactis_AD011, 'FRUt2r',...
'reactionFormula', 'fru[e] + h[e]   <=>   fru[c] + h[c]');

% Clostridium M62_1:
Clostridium_sp_M62_1 = addReaction(Clostridium_sp_M62_1, 'EX_fru(e)',...
'reactionFormula', 'fru[e]   <=>');
Clostridium_sp_M62_1 = addReaction(Clostridium_sp_M62_1, 'FRUpts',...
'reactionFormula', 'fru[e] + pep[c]   ->   f1p[c] + pyr[c]');


%% Paracasei
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'EX_kestottr(e)',...
'reactionFormula', 'kestottr[e] <=>');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'KESTOTTRASEe',...
'reactionFormula', '3.0 h2o[e] + kestottr[e] -> 3.0 fru[e] + glc_D[e]');

%% innocuum HGF2
Clostridium_innocuum_2959 = addReaction(Clostridium_innocuum_2959, 'EX_kestottr(e)',...
'reactionFormula', 'kestottr[e] <=>');
Clostridium_innocuum_2959 = addReaction(Clostridium_innocuum_2959, 'KESTOTTRASEe',...
'reactionFormula', '3.0 h2o[e] + kestottr[e] -> 3.0 fru[e] + glc_D[e]');

Clostridium_innocuum_2959 = addReaction(Clostridium_innocuum_2959, 'EX_kestopt(e)',...
'reactionFormula', 'kestopt[e] <=>');
Clostridium_innocuum_2959 = addReaction(Clostridium_innocuum_2959, 'KESTOPTASEe',...
'reactionFormula', '4.0 h2o[e] + kestopt[e] -> 4.0 fru[e] + glc_D[e]');

%% Reaction per bacterium, according to bioinformatic analysis
% L. paracasei M38
% Exoinulinase 3.2.1.80 activity (also has some invertase activity):
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'EX_inulin(e)',...
'reactionFormula', 'inulin[e] <=>');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'INULINabc',...
'reactionFormula', 'atp[c] + h2o[c] + inulin[e] -> adp[c] + h[c] + inulin[c] + pi[c]');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'INULINASE',...
'reactionFormula', '29.0 h2o[c] + inulin[c] -> 29.0 fru[c] + glc_D[c]');
% Activity Invertase 3.2.1.26 (assuming it can degrade small-chain FOS):
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'EX_kestopt(e)',...
'reactionFormula', 'kestopt[e] <=>');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = addReaction(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'KESTOPTASEe',...
'reactionFormula', '4.0 h2o[e] + kestopt[e] -> 4.0 fru[e] + glc_D[e]');

