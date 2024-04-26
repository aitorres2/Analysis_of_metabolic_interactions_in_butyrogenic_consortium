%% %% SteadyCom pipeline
clearvars
clc
%% Models load AGORA or VMH only
% initCobraToolbox(false)
changeCobraSolver('gurobi')

Coprococcus_eutactus_ATCC_27759 = readCbModel('Coprococcus_eutactus_ATCC_27759');
Bifidobacterium_animalis_subsp_lactis_PT33 =readCbModel('Bifidobacterium_animalis_subsp_lactis_PT33');
Bacteroides_thetaiotaomicron_VPI_5482 =readCbModel('Bacteroides_thetaiotaomicron_VPI_5482');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 =readCbModel('Lactobacillus_paracasei_subsp_paracasei_ATCC_25302');

models = {Coprococcus_eutactus_ATCC_27759; Bifidobacterium_animalis_subsp_lactis_PT33; Bacteroides_thetaiotaomicron_VPI_5482; Lactobacillus_paracasei_subsp_paracasei_ATCC_25302};
Names = {'Coprococcus_eutactus_ATCC_27759'; 'Bifidobacterium_animalis_subsp_lactis_PT33'; 'Bacteroides_thetaiotaomicron_VPI_5482'; 'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302'};

%% Forcing that models are not missing any necessary reaction to grow
% This function should only be used in VMH models that are not able to grow in the used medium
% Therefore, this function was only used with the Clostridium M62_1 and Clostridium
% innocuum 2959 VMH models
% 
% cultureMedium = 'cultureMediumZMBinulin_NOinulin.xlsx';
% carbonSource = 'EX_inulin(e)';
% 
% [Coprococcus_eutactus_ATCC_27759, Coprococcus_eutactus_ATCC_27759_ExcsAdded] = ...
%     loadCultureMedium(Coprococcus_eutactus_ATCC_27759, cultureMedium, carbonSource);
% 
% [Bifidobacterium_animalis_subsp_lactis_PT33, Bifidobacterium_animalis_subsp_lactis_PT33_ExcsAdded] = ...
%     loadCultureMedium(Bifidobacterium_animalis_subsp_lactis_PT33, cultureMedium, carbonSource);
% 
% [Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482_ExcsAdded] = ...
%     loadCultureMedium(Bacteroides_thetaiotaomicron_VPI_5482, cultureMedium, carbonSource);
% 
% [Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_ExcsAdded] = ...
%     loadCultureMedium(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, cultureMedium, carbonSource);

%% Condition of anaerobiosis
Coprococcus_eutactus_ATCC_27759=changeRxnBounds(Coprococcus_eutactus_ATCC_27759,'EX_o2(e)',0,'l');
Bifidobacterium_animalis_subsp_lactis_PT33=changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33,'EX_o2(e)',0,'l');
Bacteroides_thetaiotaomicron_VPI_5482=changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482,'EX_o2(e)',0,'l');
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302=changeRxnBounds(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302,'EX_o2(e)',0,'l');

models = {Coprococcus_eutactus_ATCC_27759; Bifidobacterium_animalis_subsp_lactis_PT33; Bacteroides_thetaiotaomicron_VPI_5482; Lactobacillus_paracasei_subsp_paracasei_ATCC_25302};

%% Definition of culture medium 
% maximum metabolite uptake rates
medioPUC = {'EX_inulin(e)';... % inulin
    'EX_kestopt(e)';... %FOS 4FRU-1GLC
    'EX_kesto(e)';... %FOS 2FRU-1GLC
    'EX_kestottr(e)';... %FOS 3FRU-1GLC
    'EX_sucr(e)';... %sucrose
    'EX_fru(e)';... %fructose
    'EX_glc_D(e)';... %glucose
    'EX_ala_L(e)';...
    'EX_arg_L(e)';...
    'EX_asn_L(e)';...
    'EX_asp_L(e)';...
    'EX_cys_L(e)';...
    'EX_gln_L(e)';...
    'EX_glu_L(e)';...
    'EX_gly(e)';...
    'EX_his_L(e)';...
    'EX_ile_L(e)';...
    'EX_leu_L(e)';....
    'EX_lys_L(e)';...
    'EX_met_L(e)';...
    'EX_phe_L(e)';...
    'EX_pro_L(e)';...
    'EX_ser_L(e)';...
    'EX_thr_L(e)';...
    'EX_trp_L(e)';...
    'EX_tyr_L(e)';...
    'EX_val_L(e)';...
    'EX_pi(e)';...
    'EX_k(e)';...
    'EX_mg2(e)';...
    'EX_fe3(e)';...
    'EX_fe2(e)';...
    'EX_na1(e)';...
    'EX_zn2(e)';...
    'EX_mn2(e)';...
    'EX_mobd(e)';...
    'EX_nh4(e)';...
    'EX_ca2(e)';...
    'EX_so4(e)';...
    'EX_cobalt2(e)';...
    'EX_cu2(e)';...
    'EX_cl(e)';...
    'EX_pnto_R(e)';...
    'EX_fol(e)';...
    'EX_nac(e)';...
    'EX_pydx(e)';...
    'EX_4abz(e)';...
    'EX_xan(e)';...
    'EX_inost(e)';...
    'EX_gthox(e)';...
    'EX_gthrd(e)';...
    'EX_btn(e)';...
    'EX_ribflv(e)';...
    'EX_thm(e)';...
    'EX_lipoate(e)';...
    'EX_pheme(e)';...
    'EX_sheme(e)';...
    'EX_ade(e)';...
    'EX_gua(e)';...
    'EX_ura(e)';...
    'EX_ddca(e)';...
    'EX_h2s(e)';...
    'EX_hxan(e)';...
    'EX_ni2(e)';...
    'EX_ocdca(e)';... 0.169280000000000
    'EX_thymd(e)';...
    'EX_spmd(e)';...
    'EX_q8(e)';...
    'EX_mqn7(e)';...
    'EX_mqn8(e)';...
    'EX_cgly(e)';...
    'EX_csn(e)';...
    'EX_thym(e)';...
    'EX_2obut(e)';...
    'EX_26dap_M(e)';... % THIS METABOLITE WAS NEEDED TO GROW INOCUM
    'EX_nmn(e)';... % THIS METABOLITE WAS NEEDED TO GROW INOCUM
    'EX_ins(e)'};

ratesPUC = [1;0.1;0.1;0.1;0.05;0.05;0.05;
    0.251761678;0.290359287;0.225000000000000;0.571946391;...
    0.035088736;0.820199195;0.820199195;0.157022092;0.200883011;...
    0.392993839;0.669317632;0.571069172;0.206146322;0.358782322;...
    0.758793908;0.445626943;0.342992391;0.092107931;0.163162621;...
    0.483347333;
    1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;0.169280000000000; %aca hay 37 mets mas
    1;1;1;1;1;1;1;1;1;1;1;1];

iZMB1 = [medioPUC; {'EX_h2o(e)';'EX_h2(e)';'EX_h(e)'}]; 
ratesiZMB1 = [ratesPUC; 1000; 1; 1];

%% Growth test
% here the models are constrained with the medium, and then the FBAs are made.
crecimiento = zeros(length(models),4);

% Change to add a different culture medium
media = {iZMB1};
    media_rates = {ratesiZMB1};
    sol_mono = {};

for i = 1:length(models)
    for j = 1:length(media)
        medium = media{j};
        rates = media_rates{j};

        model_test = changeRxnBounds(models{i}, models{i}.rxns(cellfun(@isempty, strfind(models{i}.rxns, 'EX_'))==0), 0, 'l');
        model_test = changeRxnBounds(model_test, medium, -rates, 'l');
        fba_model= optimizeCbModel(model_test);
        exc_idx = cellfun(@isempty, strfind(model_test.rxns, 'EX_'))==0;
        
        Coprococcus_eutactus_ATCC_27759 = changeRxnBounds(Coprococcus_eutactus_ATCC_27759, Coprococcus_eutactus_ATCC_27759.rxns(cellfun(@isempty, strfind(Coprococcus_eutactus_ATCC_27759.rxns, 'EX_'))==0), 0, 'l');
        Coprococcus_eutactus_ATCC_27759 = changeRxnBounds(Coprococcus_eutactus_ATCC_27759, medium, -rates, 'l');    
        fba_model_Coprococcus_eutactus_ATCC_27759= optimizeCbModel(Coprococcus_eutactus_ATCC_27759, 'max');
        exc_idx_Coprococcus_eutactus_ATCC_27759 = cellfun(@isempty, strfind(Coprococcus_eutactus_ATCC_27759.rxns, 'EX_'))==0;
        
        Bifidobacterium_animalis_subsp_lactis_PT33 = changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33, Bifidobacterium_animalis_subsp_lactis_PT33.rxns(cellfun(@isempty, strfind(Bifidobacterium_animalis_subsp_lactis_PT33.rxns, 'EX_'))==0), 0, 'l');
        Bifidobacterium_animalis_subsp_lactis_PT33 = changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33, medium, -rates, 'l');
        fba_model_Bifidobacterium_animalis_subsp_lactis_PT33= optimizeCbModel(Bifidobacterium_animalis_subsp_lactis_PT33, 'max');
        exc_idx_Bifidobacterium_animalis_subsp_lactis_PT33 = cellfun(@isempty, strfind(Bifidobacterium_animalis_subsp_lactis_PT33.rxns, 'EX_'))==0; 
        
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482.rxns(cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0), 0, 'l');
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, medium, -rates, 'l');    
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482= optimizeCbModel(Bacteroides_thetaiotaomicron_VPI_5482, 'max');
        exc_idx_Bacteroides_thetaiotaomicron_VPI_5482 = cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0;
                
        Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = changeRxnBounds(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.rxns(cellfun(@isempty, strfind(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.rxns, 'EX_'))==0), 0, 'l');
        Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = changeRxnBounds(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, medium, -rates, 'l');    
        fba_model_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302= optimizeCbModel(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, 'max');
        exc_idx_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = cellfun(@isempty, strfind(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.rxns, 'EX_'))==0;
               
        Coprococcus_eutactus_ATCC_27759_mono_table = table(Coprococcus_eutactus_ATCC_27759.rxns(find(fba_model_Coprococcus_eutactus_ATCC_27759.v(:,1) ~= 0),1), ...
        fba_model_Coprococcus_eutactus_ATCC_27759.v(find(fba_model_Coprococcus_eutactus_ATCC_27759.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});      
        Bifidobacterium_animalis_subsp_lactis_PT33_mono_table = table(Bifidobacterium_animalis_subsp_lactis_PT33.rxns(find(fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(:,1) ~= 0),1), ...
        fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(find(fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});    
        Bacteroides_thetaiotaomicron_VPI_5482_mono_table = table(Bacteroides_thetaiotaomicron_VPI_5482.rxns(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ...
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
        Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_mono_table = table(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.rxns(find(fba_model_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.v(:,1) ~= 0),1), ...
        fba_model_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.v(find(fba_model_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
    
        writetable(Coprococcus_eutactus_ATCC_27759_mono_table,'Coprococcus_eutactus_ATCC_27759_monocultivo_agora.xlsx');
        writetable(Bifidobacterium_animalis_subsp_lactis_PT33_mono_table,'Bifidobacterium_animalis_subsp_lactis_PT33_monocultivo_agora.xlsx');
        writetable(Bacteroides_thetaiotaomicron_VPI_5482_mono_table,'Bacteroides_thetaiotaomicron_VPI_5482_monocultivo_agora.xlsx');
        writetable(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_mono_table,'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_monocultivo_agora.xlsx');

        crecimiento(i, j) = fba_model.f;
        sol_mono{end+1} = fba_model.x;
    end   
end

%% Graphical fluxes of external reactions/metabolites of interest in monocultures
targets = {'EX_inulin(e)'; 'EX_kesto(e)'; 'EX_kestopt(e)'; 'EX_kestottr(e)'; 'EX_biomass(e)'; 'EX_ac(e)'; 'EX_but(e)'; 'EX_for(e)'; 'EX_lac_D(e)'; 'EX_lac_L(e)'; 'EX_succ(e)'; 'EX_ppa(e)'};

ids_Coprococcus_eutactus_ATCC_27759 = findRxnIDs(Coprococcus_eutactus_ATCC_27759, targets);
targets_ids_Coprococcus_eutactus_ATCC_27759 = ids_Coprococcus_eutactus_ATCC_27759(ids_Coprococcus_eutactus_ATCC_27759 > 0);
targets_filt_Coprococcus_eutactus_ATCC_27759 = targets(ids_Coprococcus_eutactus_ATCC_27759 > 0);

ids_Bifidobacterium_animalis_subsp_lactis_PT33 = findRxnIDs(Bifidobacterium_animalis_subsp_lactis_PT33, targets);
targets_ids_Bifidobacterium_animalis_subsp_lactis_PT33 = ids_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);
targets_filt_Bifidobacterium_animalis_subsp_lactis_PT33 = targets(ids_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);

ids_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(Bacteroides_thetaiotaomicron_VPI_5482, targets);
targets_ids_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Bacteroides_thetaiotaomicron_VPI_5482(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Bacteroides_thetaiotaomicron_VPI_5482 = targets(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = findRxnIDs(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, targets);
targets_ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302(ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 > 0);
targets_filt_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = targets(ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 > 0);

metabolitos_mono_Coprococcus_eutactus_ATCC_27759_table = table(targets_filt_Coprococcus_eutactus_ATCC_27759, ...
    fba_model_Coprococcus_eutactus_ATCC_27759.x(targets_ids_Coprococcus_eutactus_ATCC_27759), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bifidobacterium_animalis_subsp_lactis_PT33_table = table(targets_filt_Bifidobacterium_animalis_subsp_lactis_PT33, ...
    fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.x(targets_ids_Bifidobacterium_animalis_subsp_lactis_PT33), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482, ...
    fba_model_Bacteroides_thetaiotaomicron_VPI_5482.x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table = table(targets_filt_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, ...
    fba_model_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302.x(targets_ids_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302), ... 
    'VariableNames',{'Metabolite','v'});

%% Creation of the steadycom community model
EcCom = createMultipleSpeciesModel(models, Names);
EcCom.csense = char('E' * ones(1,numel(EcCom.mets))); 
[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom,Names);

%% Add biomass reaction information
rxnBiomass = EcCom.rxns(find(contains(EcCom.rxns, strcat(Names, 'biomass'))))
rxnBiomassId = findRxnIDs(EcCom, rxnBiomass);
EcCom.infoCom.spBm = rxnBiomass; 
EcCom.indCom.spBm = rxnBiomassId;

%% Definition of model additional constraint
% The steps to follow are:
% 1. Add the limits for the reactions of the bacteria separately
% 2. Add the limits for the reactions of the extracellular medium
% We will use iZMB1 to optimize

iZMB2 = strrep(iZMB1, '(e)', '[u]');
% Reaction restrictions for each bacterium
for i=1:length(Names)
% We set the flows of each bacterium in an unrestricted way so as not to limit the exchange between the pairs
    EcCom = changeRxnBounds(EcCom, EcCom.infoCom.EXsp(:,i),-100, 'l');
    
end

% Here we set the bounds of our culture medium in the extracellular space
EcCom = changeRxnBounds(EcCom, EcCom.rxns(startsWith(EcCom.rxns, 'EX_')), 0, 'l');
EcCom = changeRxnBounds(EcCom, iZMB2, -ratesiZMB1, 'l');
EcCom = changeRxnBounds(EcCom, EcCom.infoCom.spBm{1}, 0.1*crecimiento(1,1), 'l')
printFluxBounds(EcCom, EcCom.infoCom.spBm(1))

%% Optimización
options = struct();
options.GRguess = 0.1; % Initial value for growth rate
options.GRtol = 1e-4; % Tolerance for growth rate 
options.algorithm = 1; %  default
tic
[sol, result] = SteadyCom(EcCom, options);
toc

Coprococcus_eutactus_ATCC_27759_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1)), ... 
    'VariableNames',{'rxn','v'});

Bifidobacterium_animalis_subsp_lactis_PT33_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2)), ... 
    'VariableNames',{'rxn','v'});

Bacteroides_thetaiotaomicron_VPI_5482_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3)), ... 
    'VariableNames',{'rxn','v'});

Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,4) ~= 0),4), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,4) ~= 0),4)), ... 
    'VariableNames',{'rxn','v'});

com_table = table(EcCom.rxns(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    sol.full(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    'VariableNames',{'rxn','v'});

%% export excel tables with consumed/produced metabolites with steadycomFBA
Coprococcus_eutactus_ATCC_27759_table(~Coprococcus_eutactus_ATCC_27759_table.v,:) = []; % https://www.mathworks.com/matlabcentral/answers/407401-delete-zero-rows-from-a-table
Bifidobacterium_animalis_subsp_lactis_PT33_table(~Bifidobacterium_animalis_subsp_lactis_PT33_table.v,:) = [];
Bacteroides_thetaiotaomicron_VPI_5482_table(~Bacteroides_thetaiotaomicron_VPI_5482_table.v,:) = [];
Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table(~Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table.v,:) = [];
com_table(~com_table.v,:) = [];

writetable(Coprococcus_eutactus_ATCC_27759_table,'Coprococcus_eutactus_ATCC_27759_consorcio_con_Bifidobacterium_animalis_subsp_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_agora.xlsx');
writetable(Bifidobacterium_animalis_subsp_lactis_PT33_table,'Bifidobacterium_animalis_subsp_lactis_PT33_consorcio_con_Coprococcus_eutactus_ATCC_27759_Bacteroides_thetaiotaomicron_VPI_5482_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_agora.xlsx');
writetable(Bacteroides_thetaiotaomicron_VPI_5482_table,'Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Coprococcus_eutactus_ATCC_27759_Bifidobacterium_animalis_subsp_lactis_PT33_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_agora.xlsx');
writetable(Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table,'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_consorcio_con_Coprococcus_eutactus_ATCC_27759_Bifidobacterium_animalis_subsp_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_agora.xlsx');
writetable(com_table,'EX_comunidad_Coprococcus_eutactus_ATCC_27759_Bifidobacterium_animalis_subsp_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_agora.xlsx');

%% Graph for community external metabolites such as butyrate
targets_com_Coprococcus_eutactus_ATCC_27759 = {'Coprococcus_eutactus_ATCC_27759IEX_inulin[u]tr';
    'Coprococcus_eutactus_ATCC_27759IEX_kesto[u]tr';
    'Coprococcus_eutactus_ATCC_27759IEX_kestopt[u]tr';
    'Coprococcus_eutactus_ATCC_27759IEX_kestottr[u]tr';
    'Coprococcus_eutactus_ATCC_27759IEX_biomass[c]tr';
    'Coprococcus_eutactus_ATCC_27759IEX_ac[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_but[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_for[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_lac_D[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_lac_L[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_succ[u]tr'; 
    'Coprococcus_eutactus_ATCC_27759IEX_ppa[u]tr'}; 
    
targets_com_Bifidobacterium_animalis_subsp_lactis_PT33 = {'Bifidobacterium_animalis_subsp_lactis_PT33IEX_inulin[u]tr';
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_kesto[u]tr';
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_kestopt[u]tr';
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_kestottr[u]tr';
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_biomass[c]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_ac[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_but[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_for[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_lac_D[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_lac_L[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_succ[u]tr'; 
    'Bifidobacterium_animalis_subsp_lactis_PT33IEX_ppa[u]tr'};

targets_com_Bacteroides_thetaiotaomicron_VPI_5482 = {'Bacteroides_thetaiotaomicron_VPI_5482IEX_inulin[u]tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_kesto[u]tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_kestopt[u]tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_kestottr[u]tr';
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_biomass[c]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_ac[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_but[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_for[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_lac_D[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_lac_L[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_succ[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_ppa[u]tr'};

targets_com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = {'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_inulin[u]tr';
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_kesto[u]tr';
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_kestopt[u]tr';
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_kestottr[u]tr';
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_biomass[c]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_ac[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_but[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_for[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_lac_D[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_lac_L[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_succ[u]tr'; 
    'Lactobacillus_paracasei_subsp_paracasei_ATCC_25302IEX_ppa[u]tr'};

targets_com = {'EX_inulin[u]';
    'EX_kesto[u]';
    'EX_kestopt[u]';
    'EX_kestottr[u]';
    'EX_biomass[c]';
    'EX_ac[u]'; 
    'EX_but[u]'; 
    'EX_for[u]'; 
    'EX_lac_D[u]'; 
    'EX_lac_L[u]'; 
    'EX_succ[u]'; 
    'EX_ppa[u]'};

ids_Com_Coprococcus_eutactus_ATCC_27759 = findRxnIDs(EcCom, targets_com_Coprococcus_eutactus_ATCC_27759);
targets_ids_Com_Coprococcus_eutactus_ATCC_27759 = ids_Com_Coprococcus_eutactus_ATCC_27759(ids_Com_Coprococcus_eutactus_ATCC_27759 > 0);
targets_filt_Com_Coprococcus_eutactus_ATCC_27759 = targets_com_Coprococcus_eutactus_ATCC_27759(ids_Com_Coprococcus_eutactus_ATCC_27759 > 0);

ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = findRxnIDs(EcCom, targets_com_Bifidobacterium_animalis_subsp_lactis_PT33);
targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);
targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = targets_com_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);

ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(EcCom, targets_com_Bacteroides_thetaiotaomicron_VPI_5482);
targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482 = targets_com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = findRxnIDs(EcCom, targets_com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302);
targets_ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302(ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 > 0);
targets_filt_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 = targets_com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302(ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302 > 0);

ids_Com = findRxnIDs(EcCom, targets_com);
targets_ids_Com = ids_Com(ids_Com > 0);
targets_filt_Com = targets_com(ids_Com > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Com_Coprococcus_eutactus_ATCC_27759, '_', '\_')), sol.full(targets_ids_Com_Coprococcus_eutactus_ATCC_27759));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')});
title('mZMB + Inulin Tetra-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33, '_', '\_')), sol.full(targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')});
title('mZMB + Inulin Tetra-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{3}, '_', '\_')});
title('mZMB + Inulin Tetra-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, '_', '\_')), sol.full(targets_ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{4}, '_', '\_')});
title('mZMB + Inulin Tetra-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com, '_', '\_')), sol.full(targets_ids_Com));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_') strrep(Names{2}, '_', '\_') strrep(Names{3}, '_', '\_') strrep(Names{4}, '_', '\_')});
title('mZMB + Inulin Tetra-culture')
hold off

metabolitos_com_Coprococcus_eutactus_ATCC_27759_table = table(targets_filt_Com_Coprococcus_eutactus_ATCC_27759, ...
    sol.full(targets_ids_Com_Coprococcus_eutactus_ATCC_27759), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bifidobacterium_animalis_subsp_lactis_PT33_table = table(targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33, ...
    sol.full(targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, ...
    sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_table = table(targets_filt_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302, ...
    sol.full(targets_ids_Com_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_comunidad_totales_table = table(targets_filt_Com, ...
    sol.full(targets_ids_Com), ... 
    'VariableNames',{'Metabolite','v'});

%% SteadyCom-FVA
% maximum percentage of total biomass 
options.optBMpercent = 100;
options.rxnNameList = strcat('X_', Names);
options.optGRpercent = [0:1:100];
options.GRtol = 1e-6;
options.threads = 0;
[fvaComMin,fvaComMax] = SteadyComFVA(EcCom, options);

%% Envelope plot for relative abundance by community growth rate
grComV = result.GRmax * options.optGRpercent / 100;
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{2}, '_', '\_');strrep(Names{3}, '_', '\_');strrep(Names{4}, '_', '\_')};
f = figure;
x_lower = 10; 
x_upper = 0; 
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(Names)
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance'); 
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');
hold off

%% FVA optimizing for metabolite of interest (butyrate)
options.optBMpercent = 100;
options.rxnNameList = EcCom.rxns(find(contains(EcCom.rxns, 'IEX_but[u]')));
options.optGRpercent = [0:1:100];
options.GRtol = 1e-6;
options.threads = 0;
[fvaComMin_B,fvaComMax_B] = SteadyComFVA(EcCom, options);

%% 
grComV = result.GRmax * options.optGRpercent / 100;
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{2}, '_', '\_')};
f = figure;
x_lower = 10; 
x_upper = 0; 
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(options.rxnNameList)
    y = [fvaComMin_B(j, :), fliplr(fvaComMax_B(j, :))]; 
    y = y(~isnan(y));
    if crecimiento(j,1) == max(crecimiento(:,1))
        if min(y) < x_lower 
            x_lower = min(y);
        end 
        if max(y) > x_upper
            x_upper = max(y);
        end 
        if max(y) > y_max
            y_max = max(y)
        end 
    end 
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

tl(1) = title('\underline{SteadyComFVA butyrate}', 'Interpreter', 'latex');
ax(1).XTick = 0:0.05:y_max;
ax(1).YTick = 0:20:400;
ax(1).YTick = 0:15:200;
lg = legend(lgLabel);
% % % lg.Box = 'off';
yl(1) = ylabel('Butyrate flux (mmol gDW-1 h-1)');
xl(1) = xlabel('Community growth rate (h^{-1})');
hold off


%% Combined figures
lgLabel1 = {strrep(Names{1}, '_', '\_'); strrep(Names{2}, '_', '\_'); strrep(Names{3}, '_', '\_'); strrep(Names{4}, '_', '\_')}; % poner todos los modelos mutantes

close all

% Create the main figure
figure('Position', [100, 100, 1200, 600])

% Define subplot positions
subplotPositions = [0.1 0.1 0.4 0.8; 0.55 0.1 0.4 0.8];

% Plot the first subplot
subplot('Position', subplotPositions(1, :))
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(Names) % de 1 a n�mero de modelos mutantes
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel1)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance'); % Relative abundance Bacteria growth rate (h^{-1})
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');

% Plot the second subplot
subplot('Position', subplotPositions(2, :))
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(options.rxnNameList) % si son dos bacterias j = 1:length(options.rxnNameList); si es 1 bacteria: j = 1:length(options.rxnNameList(1))
    y = [fvaComMin_B(j, :), fliplr(fvaComMax_B(j, :))];
    y_max = max(y);
    y = y(~isnan(y));
    if crecimiento(j, 1) == max(crecimiento(:, 1))
        if min(y) < x_lower
            x_lower = min(y);
        end
        if max(y) >= x_upper
            x_upper = y_max;
        end
        if max(y) >= y_max
            y_max = max(y);
        end
    end
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

tl(1) = title('\underline{SteadyComFVA butyrate}', 'Interpreter', 'latex');
lgLabel2 = {strrep(Names{1}, '_', '\_')}; % poner a la bacteria productora de butirato, nombre2
lg = legend(lgLabel2);
yl(1) = ylabel('Butyrate flux (mmol gDW^{-1} h^{-1})');
xl(1) = xlabel('Community growth rate (h^{-1})');

% Save the combined figure
saveas(p(j, 1),'EX_combinado_Coprococcus_eutactus_ATCC_27759_Bifidobacterium_animalis_subsp_lactis_PT33_Bacteroides_thetaiotaomicron_VPI_5482_Lactobacillus_paracasei_subsp_paracasei_ATCC_25302_agora.png');
hold off
