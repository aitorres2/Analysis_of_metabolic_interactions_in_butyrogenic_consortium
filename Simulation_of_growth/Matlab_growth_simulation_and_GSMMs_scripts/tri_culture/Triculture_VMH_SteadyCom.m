%% Para "automatizar" el script hacer click en ctrl+f y reemplazar Clostridium_innocuum_HFG2 por el modelo que se desee, y lo mismo para Bacteroides_thetaiotaomicron_VPI_5482
clearvars
clc

%% Models load AGORA or VMH only
% initCobraToolbox(false)
changeCobraSolver('gurobi')

Clostridium_innocuum_HFG2 = readCbModel('Clostridium_innocuum_HFG2');
Bacteroides_thetaiotaomicron_VPI_5482 =readCbModel('Bacteroides_thetaiotaomicron_VPI_5482');
Bifidobacterium_animalis_subsp_lactis_PT33 =readCbModel('Bifidobacterium_animalis_subsp_lactis_PT33');

models = {Clostridium_innocuum_HFG2; Bacteroides_thetaiotaomicron_VPI_5482; Bifidobacterium_animalis_subsp_lactis_PT33};
Names = {'Clostridium_innocuum_HFG2'; 'Bacteroides_thetaiotaomicron_VPI_5482'; 'Bifidobacterium_animalis_subsp_lactis_PT33'};

%% Forcing that models are not missing any necessary reaction to grow
% This function should only be used in VMH models that are not able to grow in the used medium
% Therefore, this function was only used with the Clostridium M62_1 and Clostridium
% innocuum 2959 VMH models
% 
% cultureMedium = 'cultureMediumZMBinulin_NOinulin.xlsx';
% carbonSource = 'EX_inulin(e)';
% 
% [Clostridium_innocuum_HFG2, Clostridium_innocuum_HFG2_ExcsAdded] = ...
%     loadCultureMedium(Clostridium_innocuum_HFG2, cultureMedium, carbonSource);
% 
% [Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482_ExcsAdded] = ...
%     loadCultureMedium(Bacteroides_thetaiotaomicron_VPI_5482, cultureMedium, carbonSource);
% 
% [Bifidobacterium_animalis_subsp_lactis_PT33, Bifidobacterium_animalis_subsp_lactis_PT33_ExcsAdded] = ...
%     loadCultureMedium(Bifidobacterium_animalis_subsp_lactis_PT33, cultureMedium, carbonSource);

%% Condition of anaerobiosis
Clostridium_innocuum_HFG2=changeRxnBounds(Clostridium_innocuum_HFG2,'EX_o2(e)',0,'l');
Bacteroides_thetaiotaomicron_VPI_5482=changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482,'EX_o2(e)',0,'l');
Bifidobacterium_animalis_subsp_lactis_PT33=changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33,'EX_o2(e)',0,'l');
models = {Clostridium_innocuum_HFG2; Bacteroides_thetaiotaomicron_VPI_5482; Bifidobacterium_animalis_subsp_lactis_PT33};

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
crecimiento = zeros(length(models),3);

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
        
        Clostridium_innocuum_HFG2 = changeRxnBounds(Clostridium_innocuum_HFG2, Clostridium_innocuum_HFG2.rxns(cellfun(@isempty, strfind(Clostridium_innocuum_HFG2.rxns, 'EX_'))==0), 0, 'l');
        Clostridium_innocuum_HFG2 = changeRxnBounds(Clostridium_innocuum_HFG2, medium, -rates, 'l');    
        fba_model_Clostridium_innocuum_HFG2= optimizeCbModel(Clostridium_innocuum_HFG2, 'max');
        exc_idx_Clostridium_innocuum_HFG2 = cellfun(@isempty, strfind(Clostridium_innocuum_HFG2.rxns, 'EX_'))==0;
        
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482.rxns(cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0), 0, 'l');
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, medium, -rates, 'l');
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482= optimizeCbModel(Bacteroides_thetaiotaomicron_VPI_5482, 'max');
        exc_idx_Bacteroides_thetaiotaomicron_VPI_5482 = cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0; 
        
        Bifidobacterium_animalis_subsp_lactis_PT33 = changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33, Bifidobacterium_animalis_subsp_lactis_PT33.rxns(cellfun(@isempty, strfind(Bifidobacterium_animalis_subsp_lactis_PT33.rxns, 'EX_'))==0), 0, 'l');
        Bifidobacterium_animalis_subsp_lactis_PT33 = changeRxnBounds(Bifidobacterium_animalis_subsp_lactis_PT33, medium, -rates, 'l');    
        fba_model_Bifidobacterium_animalis_subsp_lactis_PT33= optimizeCbModel(Bifidobacterium_animalis_subsp_lactis_PT33, 'max');
        exc_idx_Bifidobacterium_animalis_subsp_lactis_PT33 = cellfun(@isempty, strfind(Bifidobacterium_animalis_subsp_lactis_PT33.rxns, 'EX_'))==0;
                
        Clostridium_innocuum_HFG2_mono_table = table(Clostridium_innocuum_HFG2.rxns(find(fba_model_Clostridium_innocuum_HFG2.v(:,1) ~= 0),1), ...
        fba_model_Clostridium_innocuum_HFG2.v(find(fba_model_Clostridium_innocuum_HFG2.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});      
        Bacteroides_thetaiotaomicron_VPI_5482_mono_table = table(Bacteroides_thetaiotaomicron_VPI_5482.rxns(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ...
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});    
        Bifidobacterium_animalis_subsp_lactis_PT33_mono_table = table(Bifidobacterium_animalis_subsp_lactis_PT33.rxns(find(fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(:,1) ~= 0),1), ...
        fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(find(fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
        
        writetable(Clostridium_innocuum_HFG2_mono_table,'Clostridium_innocuum_HFG2_monocultivo_agora.xlsx');
        writetable(Bacteroides_thetaiotaomicron_VPI_5482_mono_table,'Bacteroides_thetaiotaomicron_VPI_5482_monocultivo_agora.xlsx');
        writetable(Bifidobacterium_animalis_subsp_lactis_PT33_mono_table,'Bifidobacterium_animalis_subsp_lactis_PT33_monocultivo_agora.xlsx');

        crecimiento(i, j) = fba_model.f;
        sol_mono{end+1} = fba_model.x;
    end   
end

%% Graphical fluxes of external reactions/metabolites of interest in monocultures
targets = {'EX_inulin(e)'; 'EX_kesto(e)'; 'EX_kestopt(e)'; 'EX_kestottr(e)'; 'EX_biomass(e)'; 'EX_ac(e)'; 'EX_but(e)'; 'EX_for(e)'; 'EX_lac_D(e)'; 'EX_lac_L(e)'; 'EX_succ(e)'; 'EX_ppa(e)'};

ids_Clostridium_innocuum_HFG2 = findRxnIDs(Clostridium_innocuum_HFG2, targets);
targets_ids_Clostridium_innocuum_HFG2 = ids_Clostridium_innocuum_HFG2(ids_Clostridium_innocuum_HFG2 > 0);
targets_filt_Clostridium_innocuum_HFG2 = targets(ids_Clostridium_innocuum_HFG2 > 0);

ids_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(Bacteroides_thetaiotaomicron_VPI_5482, targets);
targets_ids_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Bacteroides_thetaiotaomicron_VPI_5482(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Bacteroides_thetaiotaomicron_VPI_5482 = targets(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Bifidobacterium_animalis_subsp_lactis_PT33 = findRxnIDs(Bifidobacterium_animalis_subsp_lactis_PT33, targets);
targets_ids_Bifidobacterium_animalis_subsp_lactis_PT33 = ids_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);
targets_filt_Bifidobacterium_animalis_subsp_lactis_PT33 = targets(ids_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);

metabolitos_mono_Clostridium_innocuum_HFG2_table = table(targets_filt_Clostridium_innocuum_HFG2, ...
    fba_model_Clostridium_innocuum_HFG2.x(targets_ids_Clostridium_innocuum_HFG2), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482, ...
    fba_model_Bacteroides_thetaiotaomicron_VPI_5482.x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bifidobacterium_animalis_subsp_lactis_PT33_table = table(targets_filt_Bifidobacterium_animalis_subsp_lactis_PT33, ...
    fba_model_Bifidobacterium_animalis_subsp_lactis_PT33.x(targets_ids_Bifidobacterium_animalis_subsp_lactis_PT33), ... 
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
for i=1:length(Names)
    EcCom = changeRxnBounds(EcCom, EcCom.infoCom.EXsp(:,i),-100, 'l');
    
end

EcCom = changeRxnBounds(EcCom, EcCom.rxns(startsWith(EcCom.rxns, 'EX_')), 0, 'l');
EcCom = changeRxnBounds(EcCom, iZMB2, -ratesiZMB1, 'l');
EcCom = changeRxnBounds(EcCom, EcCom.infoCom.spBm{1}, 0.1*crecimiento(1,1), 'l')
printFluxBounds(EcCom, EcCom.infoCom.spBm(1))

%% Optimization
changeCobraSolver('gurobi')
options = struct();
options.GRguess = 0.1; % Initial value for growth rate
options.GRtol = 1e-4; % Tolerance for growth rate 
options.algorithm = 1; %  default
tic
[sol, result] = SteadyCom(EcCom, options);
toc

Clostridium_innocuum_HFG2_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1)), ... 
    'VariableNames',{'rxn','v'});

Bacteroides_thetaiotaomicron_VPI_5482_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2)), ... 
    'VariableNames',{'rxn','v'});

Bifidobacterium_animalis_subsp_lactis_PT33_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,3) ~= 0),3)), ... 
    'VariableNames',{'rxn','v'});

com_table = table(EcCom.rxns(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    sol.full(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    'VariableNames',{'rxn','v'});

%% export excel tables with consumed/produced metabolites with steadycomFBA
Clostridium_innocuum_HFG2_table(~Clostridium_innocuum_HFG2_table.v,:) = []; % https://www.mathworks.com/matlabcentral/answers/407401-delete-zero-rows-from-a-table
Bacteroides_thetaiotaomicron_VPI_5482_table(~Bacteroides_thetaiotaomicron_VPI_5482_table.v,:) = [];
Bifidobacterium_animalis_subsp_lactis_PT33_table(~Bifidobacterium_animalis_subsp_lactis_PT33_table.v,:) = [];
com_table(~com_table.v,:) = [];

writetable(Clostridium_innocuum_HFG2_table,'Clostridium_innocuum_HFG2_consorcio_con_Bacteroides_thetaiotaomicron_VPI_5482_Bifidobacterium_animalis_subsp_lactis_PT33_agora.xlsx');
writetable(Bacteroides_thetaiotaomicron_VPI_5482_table,'Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Clostridium_innocuum_HFG2_Bifidobacterium_animalis_subsp_lactis_PT33_agora.xlsx');
writetable(Bifidobacterium_animalis_subsp_lactis_PT33_table,'Bifidobacterium_animalis_subsp_lactis_PT33_consorcio_con_Clostridium_innocuum_HFG2_Bacteroides_thetaiotaomicron_VPI_5482_agora.xlsx');
writetable(com_table,'EX_comunidad_Clostridium_innocuum_HFG2_Bacteroides_thetaiotaomicron_VPI_5482_Bifidobacterium_animalis_subsp_lactis_PT33_agora.xlsx');

%% Graph for community external metabolites such as butyrate
targets_com_Clostridium_innocuum_HFG2 = {'Clostridium_innocuum_HFG2IEX_inulin[u]tr';
    'Clostridium_innocuum_HFG2IEX_kesto[u]tr';
    'Clostridium_innocuum_HFG2IEX_kestopt[u]tr';
    'Clostridium_innocuum_HFG2IEX_kestottr[u]tr';
    'Clostridium_innocuum_HFG2IEX_biomass[c]tr';
    'Clostridium_innocuum_HFG2IEX_ac[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_but[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_for[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_lac_D[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_lac_L[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_succ[u]tr'; 
    'Clostridium_innocuum_HFG2IEX_ppa[u]tr'}; 
    
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

ids_Com_Clostridium_innocuum_HFG2 = findRxnIDs(EcCom, targets_com_Clostridium_innocuum_HFG2);
targets_ids_Com_Clostridium_innocuum_HFG2 = ids_Com_Clostridium_innocuum_HFG2(ids_Com_Clostridium_innocuum_HFG2 > 0);
targets_filt_Com_Clostridium_innocuum_HFG2 = targets_com_Clostridium_innocuum_HFG2(ids_Com_Clostridium_innocuum_HFG2 > 0);

ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(EcCom, targets_com_Bacteroides_thetaiotaomicron_VPI_5482);
targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482 = targets_com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = findRxnIDs(EcCom, targets_com_Bifidobacterium_animalis_subsp_lactis_PT33);
targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);
targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33 = targets_com_Bifidobacterium_animalis_subsp_lactis_PT33(ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33 > 0);

ids_Com = findRxnIDs(EcCom, targets_com);
targets_ids_Com = ids_Com(ids_Com > 0);
targets_filt_Com = targets_com(ids_Com > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Com_Clostridium_innocuum_HFG2, '_', '\_')), sol.full(targets_ids_Com_Clostridium_innocuum_HFG2));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')});
title('mZMB + Inulin Tri-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')});
title('mZMB + Inulin Tri-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33, '_', '\_')), sol.full(targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{3}, '_', '\_')});
title('mZMB + Inulin Tri-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com, '_', '\_')), sol.full(targets_ids_Com));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_') strrep(Names{2}, '_', '\_') strrep(Names{3}, '_', '\_')});
title('mZMB + Inulin Tri-culture')
hold off

metabolitos_com_Clostridium_innocuum_HFG2_table = table(targets_filt_Com_Clostridium_innocuum_HFG2, ...
    sol.full(targets_ids_Com_Clostridium_innocuum_HFG2), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, ...
    sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bifidobacterium_animalis_subsp_lactis_PT33_table = table(targets_filt_Com_Bifidobacterium_animalis_subsp_lactis_PT33, ...
    sol.full(targets_ids_Com_Bifidobacterium_animalis_subsp_lactis_PT33), ... 
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
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{2}, '_', '\_');strrep(Names{3}, '_', '\_')};
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
ax(1).YTick = 0:2:50;
lg = legend(lgLabel);
% % % lg.Box = 'off';
yl(1) = ylabel('Butyrate flux (mmol gDW-1 h-1)');
xl(1) = xlabel('Community growth rate (h^{-1})');
hold off

%% Combined figures
lgLabel1 = {strrep(Names{1}, '_', '\_'); strrep(Names{2}, '_', '\_'); strrep(Names{3}, '_', '\_')}; % poner todos los modelos mutantes

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
saveas(p(j, 1),'EX_combinado_Clostridium_innocuum_HFG2_Bacteroides_thetaiotaomicron_VPI_5482_Bifidobacterium_animalis_subsp_lactis_PT33_agora.png');
hold off
