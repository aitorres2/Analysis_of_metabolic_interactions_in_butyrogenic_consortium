%% refinement of gapseq models for simulation on inulin
% Change working directory if necessary

% initCobraToolbox(false)
changeCobraSolver('gurobi')

%% load draft models from gapseq
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

%% Establish biomass as objective function
objectiveAbbr_Bifidobacterium_animalis_lactis_PT33 = checkObjective(Bifidobacterium_animalis_lactis_PT33);
objectiveAbbr_Clostridium_sp_7_2_43FAA = checkObjective(Clostridium_sp_7_2_43FAA);
objectiveAbbr_Clostridium_sp_M62 = checkObjective(Clostridium_sp_M62);
%%
Bifidobacterium_animalis_lactis_PT33 = changeObjective(Bifidobacterium_animalis_lactis_PT33, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1
objectiveAbbr_Bifidobacterium_animalis_lactis_PT33 = checkObjective(Bifidobacterium_animalis_lactis_PT33);

%%
objectiveAbbr_Bifidobacterium_longum_PT8 = checkObjective(Bifidobacterium_longum_PT8);
objectiveAbbr_Clostridium_innocuum_HGF2_Manual_Inulin = checkObjective(Clostridium_innocuum_HGF2_Manual_Inulin);
%%
% cpd00023 = L-Glutamate
% cpd00281 = GABA
% cpd00211 = butyrate
% cpd00029 = acetate
% cpd11602 = inulin
% cpd00036 = succinate
% cpd00047 = formate
% cpd00159 = L-lactate
% cpd00221 = D-lactate
% cpd01022 = lactate
% cpd00141 = propionate
% cpd00363 = ethanol
%     'EX_cpd02298_e0';... FOS
%     'EX_cpd00027_e0';... % D-glucose
%     'EX_cpd00082_e0';... % D-fructose
%     'EX_cpd00076_e0';... % Sucrose
    
%%
Bifidobacterium_longum_PT8 = changeObjective(Bifidobacterium_longum_PT8, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1
Clostridium_innocuum_HGF2_Manual_Inulin = changeObjective(Clostridium_innocuum_HGF2_Manual_Inulin, 'bio1_biomass', 1);  % Set the biomass reaction coefficient to 1

%% THIS REACTION SHOULD NOT BE IN THE MODEL
% look for what reactions the model should not have or should have:
%%
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00211[e0]"); % butyrate; no deberia tener %%%% no tiene
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd11602[e0]"); % inulina; no deberia tener %%%% no tiene

%%
surfNet(Clostridium_sp_7_2_43FAA, "cpd00211[e0]"); % butyrate; si deberia tener %%%% si tiene
%%
surfNet(Clostridium_sp_M62, "cpd00211[e0]"); % butyrate; si deberia tener
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

%% Adding inulin and FOS modelseed reactions
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'EX_cpd11602_e0',...
'reactionFormula', 'cpd11602[e0] <=>');
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'EX_cpd02298_e0',...
'cpd02298[e0] <=>');

%Reacción de transporte de inulina con protones: rxn05695
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn05695_e0',...
'reactionFormula', 'cpd00067[c0] + cpd11602[c0] <=> cpd00067[e0] + cpd11602[e0]');

%Descomposición de 1F-beta-D-Fructosylsucrose (repetida, por lo que sólo
%listo una): rxn23732 / rxn23733
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23732_e0',...
'reactionFormula', 'cpd00076[e0] + cpd11602[e0] => 2 cpd02298[e0]');

%Transporte simple de inulina: rxn28593
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn28593_e0',...
'reactionFormula', 'cpd11602[c0] <=> cpd11602[e0]');

%Reacción de consumo de ATP con transporte de inulina: rxn29776
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn29776_e0',...
'reactionFormula', 'cpd00001[e0] + cpd00002[e0] + cpd11602[c0] => cpd00008[e0] + cpd00009[e0] + cpd00067[e0] + cpd11602[e0]');

%Reacción de intercambio de protones con transporte de inulina (en ambas
%direcciones): rxn29847
Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn29847_e0',...
'reactionFormula', 'cpd00067[e0] + cpd11602[c0] <=> cpd00067[c0] + cpd11602[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn00012_e0',...
'reactionFormula', '2 cpd00076[e0]  <=> cpd00027[e0] + cpd02298[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23739_e0',...
'reactionFormula', 'cpd00076[e0] + cpd02298[e0]  <=> cpd00190[e0] + cpd22431[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn23749_e0',...
'reactionFormula', 'cpd00001[e0] + cpd02298[e0]  <=> cpd00076[e0] + cpd19102[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn29739_e0',...
'cpd00001[e0] + cpd00002[e0] + cpd02298[c0] => cpd00008[e0] + cpd00009[e0] + cpd00067[e0] + cpd02298[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn37232_e0',...
'reactionFormula', 'cpd00001[e0] + cpd02298[e0]  -> cpd00076[e0] + cpd30321[e0]');

Clostridium_innocuum_HFG2 = addReaction(Clostridium_innocuum_HFG2, 'rxn38360_e0',...
'2 cpd00076[e0]  <=> cpd00190[e0] + cpd02298[e0]');

%% Elimnar butyrate reactions in bacteroids and paracasei
% surfNet(Clostridium_innocuum_HFG2, "cpd00211[e0]"); % inulina %%%% no tiene este metabolito
surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "cpd00211[e0]"); % inulina %%%% no tiene este metabolito
Bacteroides_thetaiotaomicron_VPI_5482 = removeRxns(Bacteroides_thetaiotaomicron_VPI_5482, 'rxn05683_c0');
Bacteroides_thetaiotaomicron_VPI_5482 = removeRxns(Bacteroides_thetaiotaomicron_VPI_5482, 'EX_cpd00211_e0');

%%
Lacticaseibacillus_paracasei_M38 = readCbModel('Lacticaseibacillus_paracasei_M38');
%%
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00211[e0]"); % butirato %%%% no tiene este metabolito
Lacticaseibacillus_paracasei_M38 = removeRxns(Lacticaseibacillus_paracasei_M38, 'rxn05683_c0');
Lacticaseibacillus_paracasei_M38 = removeRxns(Lacticaseibacillus_paracasei_M38, 'EX_cpd00211_e0');
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00211[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd11602[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00027[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00190[e0]"); % butirato %%%% no tiene este metabolito
surfNet(Lacticaseibacillus_paracasei_M38, "cpd00190[e0]"); % butirato %%%% no tiene este metabolito

%% Adding inulin and FOS modelseed reactions
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'EX_cpd11602_e0',...
'reactionFormula', 'cpd11602[e0] <=>');
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'EX_cpd02298_e0',...
'cpd02298[e0] <=>');

%Reacción de transporte de inulina con protones: rxn05695
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn05695_e0',...
'reactionFormula', 'cpd00067[c0] + cpd11602[c0] <=> cpd00067[e0] + cpd11602[e0]');

%Descomposición de 1F-beta-D-Fructosylsucrose (repetida, por lo que sólo
%listo una): rxn23732 / rxn23733
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23732_e0',...
'reactionFormula', 'cpd00076[e0] + cpd11602[e0] => 2 cpd02298[e0]');

%Transporte simple de inulina: rxn28593
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn28593_e0',...
'reactionFormula', 'cpd11602[c0] <=> cpd11602[e0]');

%Reacción de consumo de ATP con transporte de inulina: rxn29776
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn29776_e0',...
'reactionFormula', 'cpd00001[e0] + cpd00002[e0] + cpd11602[c0] => cpd00008[e0] + cpd00009[e0] + cpd00067[e0] + cpd11602[e0]');

%Reacción de intercambio de protones con transporte de inulina (en ambas
%direcciones): rxn29847
Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn29847_e0',...
'reactionFormula', 'cpd00067[e0] + cpd11602[c0] <=> cpd00067[c0] + cpd11602[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn00012_e0',...
'reactionFormula', '2 cpd00076[e0]  <=> cpd00027[e0] + cpd02298[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23739_e0',...
'reactionFormula', 'cpd00076[e0] + cpd02298[e0]  <=> cpd00190[e0] + cpd22431[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn23749_e0',...
'reactionFormula', 'cpd00001[e0] + cpd02298[e0]  <=> cpd00076[e0] + cpd19102[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn29739_e0',...
'cpd00001[e0] + cpd00002[e0] + cpd02298[c0] => cpd00008[e0] + cpd00009[e0] + cpd00067[e0] + cpd02298[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn37232_e0',...
'reactionFormula', 'cpd00001[e0] + cpd02298[e0]  -> cpd00076[e0] + cpd30321[e0]');

Lacticaseibacillus_paracasei_M38 = addReaction(Lacticaseibacillus_paracasei_M38, 'rxn38360_e0',...
'2 cpd00076[e0]  <=> cpd00190[e0] + cpd02298[e0]');


%%
% cpd00023 = L-Glutamate
% cpd00281 = GABA
% cpd00211 = butyrate
% cpd00029 = acetate
% cpd11602 = inulin
% cpd00036 = succinate
% cpd00047 = formate
% cpd00159 = L-lactate
% cpd00221 = D-lactate
% cpd01022 = lactate
% cpd00141 = propionate
% cpd00363 = ethanol
%     'EX_cpd02298_e0';... FOS
%     'EX_cpd00027_e0';... % D-glucose
%     'EX_cpd00082_e0';... % D-fructose
%     'EX_cpd00076_e0';... % Sucrose

surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00082[e0]"); % fructose %%%% does not have this metabolite
%%
surfNet(Bacteroides_thetaiotaomicron_VPI_5482, "cpd00221[e0]"); % butyrate %%%% does not have this metabolite

%%

surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00082[c0]"); % fructose %%%% does not have this metabolite

%%
surfNet(Clostridium_sp_M62, "cpd00082[e0]"); % %%%% does not have this metabolite
%%
surfNet(Bifidobacterium_animalis_lactis_PT33, "cpd00076[e0]"); % %%%% does not have this metabolite
surfNet(Clostridium_sp_M62, "cpd00076[e0]"); % %%%% does not have this metabolite

%% Seguimiento rutas de metabolitos
% https://opencobra.github.io/cobratoolbox/latest/modules/analysis/exploration/index.html?highlight=biomassPrecursorCheck#src.analysis.exploration.biomassPrecursorCheck

% biomassPrecursorCheck(Lacticaseibacillus_paracasei_M38) %la consola me imprime los metabolitos esencailes apra la biomasa pero q no los pouede producir el modelo

% cpd00023 = L-Glutamate
% cpd00281 = GABA
% cpd00211 = butyrate
% cpd00029 = acetate
% cpd11602 = inulin
% cpd00036 = succinate
% cpd00047 = formate
% cpd00159 = L-lactate
% cpd00221 = D-lactate
% cpd01022 = lactate
% cpd00141 = propionate
% cpd00363 = ethanol
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
% 
% % Comparison with AGORA models
% 
% Bifidobacterium_longum_DJO10A = readCbModel('Bifidobacterium_longum_DJO10A');
% Lactobacillus_casei_ATCC_334 =readCbModel('Lactobacillus_casei_ATCC_334');
% 
% surfNet(Bifidobacterium_longum_DJO10A, "ac[e]");% cpd00029 = acetate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00036[e0]");% cpd00036 = succinate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00047[e0]");% cpd00047 = formate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00159[e0]");% cpd00159 = L-lactate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00221[e0]");% cpd00221 = D-lactate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd01022[e0]");% cpd01022 = lactate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00141[e0]");% cpd00141 = propionate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00141[e0]");% cpd00211 = butyrate
% surfNet(Bifidobacterium_longum_DJO10A, "cpd00363[e0]");% cpd00363 = ethanol
% 
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00029[e0]");% cpd00029 = acetate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00036[e0]");% cpd00036 = succinate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00047[e0]");% cpd00047 = formate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00159[e0]");% cpd00159 = L-lactate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00221[e0]");% cpd00221 = D-lactate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd01022[e0]");% pd01022 = lactate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00141[e0]");% cpd00141 = propionate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00211[e0]");% cpd00211 = butyrate
% surfNet(Lactobacillus_casei_ATCC_334, "cpd00363[e0]");% cpd00363 = ethanol
