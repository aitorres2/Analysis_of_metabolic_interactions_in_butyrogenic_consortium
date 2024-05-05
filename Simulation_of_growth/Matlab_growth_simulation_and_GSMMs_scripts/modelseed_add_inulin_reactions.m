%% refinamiento modelos gapseq para simulacion en inulina
% /media/alexis/hdd2/objetivo_1_tesis_doctoral/matlab_scripts_genomas

% initCobraToolbox(false)

changeCobraSolver('gurobi')

%% cargar modelos recien salidos de gapseq
Bifidobacterium_animalis_lactis_PT33  = readCbModel('Bifidobacterium_animalis_lactis_PT33');
Clostridium_sp_7_2_43FAA = readCbModel('Clostridium_sp_7_2_43FAA');
Clostridium_sp_M62  = readCbModel('Clostridium_sp_M62');

%% PT33
Bifidobacterium_animalis_lactis_PT33  = readCbModel('Bifidobacterium_animalis_lactis_PT33');

%%
% Define the old reaction name and the new reaction name
oldReactionName = 'bio1';
newReactionName = 'bio1_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Bifidobacterium_animalis_lactis_PT33.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Bifidobacterium_animalis_lactis_PT33.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%%
% Define the old reaction name and the new reaction name
oldReactionName = 'bio1';
newReactionName = 'bio1_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_sp_7_2_43FAA.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_sp_7_2_43FAA.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end

%%
% Define the old reaction name and the new reaction name
oldReactionName = 'bio1';
newReactionName = 'bio1_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_sp_M62.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_sp_M62.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end

%%
% Define the old reaction name and the new reaction name
oldReactionName = 'bio1';
newReactionName = 'bio1_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Bifidobacterium_longum_PT8.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Bifidobacterium_longum_PT8.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end
%%
% Define the old reaction name and the new reaction name
oldReactionName = 'bio1';
newReactionName = 'bio1_biomass';

% Find the index of the reaction with the old name
reactionIndex = find(strcmp(Clostridium_innocuum_HGF2_Manual_Inulin.rxns, oldReactionName));

% Check if the reaction is found
if ~isempty(reactionIndex)
    % Change the name of the reaction
    Clostridium_innocuum_HGF2_Manual_Inulin.rxns{reactionIndex} = newReactionName;
    disp(['Reaction name changed from ' oldReactionName ' to ' newReactionName]);
else
    disp(['Reaction ' oldReactionName ' not found in the model.']);
end

%%
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00211[e0]"); % butyrate; no deberia tener %%%% no tiene
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd11602[e0]"); % inulina; no deberia tener %%%% no tiene
%%
surfNet(Clostridium_sp_7_2_43FAA, "cpd00211[e0]"); % butyrate; 
%%
surfNet(Clostridium_sp_M62, "cpd00211[e0]"); % butyrate; 
%%


%%
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn05695_e0')
%%
% printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23732_c0')

clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23732_c0')

%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn28593_e0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn29776_e0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn00012_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23732_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23734_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23739_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn23749_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn28566_e0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn29535_e0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn29739_e0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn37232_c0')
%%
clc
printRxnFormula(Clostridium_innocuum_HGF2_Manual_Inulin, 'rxn38360_c0')

%% Añadir reacciones de inulina y FOS modelseed

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn28593_e0',...
'reactionFormula', 'cpd11602[e0] <=>');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn28566_e0',...
'cpd02298[e0] <=>');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn29776_e0',...
'reactionFormula', 'cpd00001[c0] + cpd00002[c0] + cpd11602[e0] -> cpd00008[c0] + cpd00009[c0] + cpd00067[c0] + cpd11602[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23732_c0',...
'reactionFormula', 'cpd00076[c0] + cpd11602[c0] -> 2 cpd02298[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn00012_c0',...
'reactionFormula', '2 cpd00076[c0]  <=> cpd00027[c0] + cpd02298[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23739_c0',...
'reactionFormula', 'cpd00076[c0] + cpd02298[c0]  <=> cpd00190[c0] + cpd22431[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23749_c0',...
'reactionFormula', 'cpd00001[c0] + cpd02298[c0]  <=> cpd00076[c0] + cpd19102[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn29739_e0',...
'cpd00001[c0] + cpd00002[c0] + cpd02298[e0] => cpd00008[c0] + cpd00009[c0] + cpd00067[c0] + cpd02298[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn37232_e0',...
'cpd00001[c0] + cpd02298[c0]  -> cpd00076[c0] + cpd30321[c0]');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn38360_e0',...
'2 cpd00076[c0]  <=> cpd00190[c0] + cpd02298[c0]');

%% Elimnar reacciones butirato en bacteroides y paracasei
% surfNet(Clostridium_innocuum_HFG2, "cpd00211[e0]"); % inulina %%%% no tiene este metabolito
surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "cpd00211[e0]"); % inulina %%%% no tiene este metabolito
Bacteroides_thetaiotaomicron_VPI_5482 = removeRxns(Bacteroides_thetaiotaomicron_VPI_5482, 'rxn05683_c0');
Bacteroides_thetaiotaomicron_VPI_5482 = removeRxns(Bacteroides_thetaiotaomicron_VPI_5482, 'EX_cpd00211_e0');

%%
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00211[e0]"); % butirato %%%% no tiene este metabolito
Lacticaseibacillus_paracasei_M38 = removeRxns(Lacticaseibacillus_paracasei_M38, 'rxn05683_c0');
Lacticaseibacillus_paracasei_M38 = removeRxns(Lacticaseibacillus_paracasei_M38, 'EX_cpd00211_e0');
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00211[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd11602[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00027[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00190[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00190[e0]"); % butirato %%%% no tiene este metabolito
%%
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn28593_e0',...
'reactionFormula', 'cpd11602[e0] <=>');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn28566_e0',...
'cpd02298[e0] <=>');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn29776_e0',...
'reactionFormula', 'cpd00001[c0] + cpd00002[c0] + cpd11602[e0] -> cpd00008[c0] + cpd00009[c0] + cpd00067[c0] + cpd11602[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23732_c0',...
'reactionFormula', 'cpd00076[c0] + cpd11602[c0] -> 2 cpd02298[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn00012_c0',...
'reactionFormula', '2 cpd00076[c0]  <=> cpd00027[c0] + cpd02298[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23739_c0',...
'reactionFormula', 'cpd00076[c0] + cpd02298[c0]  <=> cpd00190[c0] + cpd22431[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23749_c0',...
'reactionFormula', 'cpd00001[c0] + cpd02298[c0]  <=> cpd00076[c0] + cpd19102[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn29739_e0',...
'cpd00001[c0] + cpd00002[c0] + cpd02298[e0] => cpd00008[c0] + cpd00009[c0] + cpd00067[c0] + cpd02298[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn37232_e0',...
'cpd00001[c0] + cpd02298[c0]  -> cpd00076[c0] + cpd30321[c0]');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn38360_e0',...
'2 cpd00076[c0]  <=> cpd00190[c0] + cpd02298[c0]');

surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00082[e0]"); % fructose %%%% no tiene este metabolito
%%
surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "cpd00221[e0]"); % butirato %%%% no tiene este metabolito

%%

surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00082[c0]"); % fructose %%%% no tiene este metabolito

%%
surfNet(Clostridium_sp_M62, "cpd00082[e0]"); % butirato %%%% no tiene este metabolito
%%
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00076[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Clostridium_sp_M62, "cpd00076[e0]"); % butirato %%%% no tiene este metabolito

% 
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00029[e0]");% cpd00029 = acetate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00036[e0]");% cpd00036 = succinate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00047[e0]");% cpd00047 = formate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00159[e0]");% cpd00159 = L-lactate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00221[e0]");% cpd00221 = D-lactate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd01022[e0]");% cpd01022 = lactate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00141[e0]");% cpd00141 = propionate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00141[e0]");% cpd00211 = butyrate
% surfNet(Lacticaseibacillus_paracasei_M38, "cpd00363[e0]");% cpd00363 = ethanol
% 
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00029[e0]");% cpd00029 = acetate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00036[e0]");% cpd00036 = succinate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00047[e0]");% cpd00047 = formate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00159[e0]");% cpd00159 = L-lactate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00221[e0]");% cpd00221 = D-lactate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd01022[e0]");% cpd01022 = lactate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00141[e0]");% cpd00141 = propionate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00211[e0]");% cpd00211 = butyrate
% surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00363[e0]");% cpd00363 = ethanol

% surfNet(Clostridium_sp_M62, "cpd00029[e0]");% cpd00029 = acetate
% surfNet(Clostridium_sp_M62, "cpd00036[e0]");% cpd00036 = succinate
% surfNet(Clostridium_sp_M62, "cpd00047[e0]");% cpd00047 = formate
% surfNet(Clostridium_sp_M62, "cpd00159[e0]");% cpd00159 = L-lactate
% surfNet(Clostridium_sp_M62, "cpd00221[e0]");% cpd00221 = D-lactate
% surfNet(Clostridium_sp_M62, "cpd01022[e0]");% pd01022 = lactate
% surfNet(Clostridium_sp_M62, "cpd00141[e0]");% cpd00141 = propionate
% surfNet(Clostridium_sp_M62, "cpd00211[e0]");% cpd00211 = butyrate
% surfNet(Clostridium_sp_M62, "cpd00363[e0]");% cpd00363 = ethanol
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00211[e0]");% cpd00211 = butyrate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00363[e0]");% cpd00363 = ethanol
