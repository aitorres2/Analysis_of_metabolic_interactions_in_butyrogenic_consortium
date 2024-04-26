%% %% SteadyCom pipeline
clearvars
clc

%% Models load AGORA or VMH only
% initCobraToolbox(false)
changeCobraSolver('gurobi')

Clostridium_symbiosum_WAL_14673 = readCbModel('Clostridium_symbiosum_WAL_14673');
Bacteroides_thetaiotaomicron_VPI_5482 = readCbModel('Bacteroides_thetaiotaomicron_VPI_5482');

%% Forcing that models are not missing any necessary reaction to grow
% This function should only be used in VMH models that are not able to grow in the used medium
% Therefore, this function was only used with the Clostridium M62_1 and Clostridium
% innocuum 2959 VMH models

% cultureMedium = 'cultureMediumZMBinulin_NOinulin.xlsx';
% carbonSource = 'EX_inulin(e)';
% 
% [Clostridium_symbiosum_WAL_14673, Clostridium_symbiosum_WAL_14673_ExcsAdded] = ...
%     loadCultureMedium(Clostridium_symbiosum_WAL_14673, cultureMedium, carbonSource);
% 
% [Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482_ExcsAdded] = ...
%     loadCultureMedium(Bacteroides_thetaiotaomicron_VPI_5482, cultureMedium, carbonSource);

%% Condition of anaerobiosis
Clostridium_symbiosum_WAL_14673=changeRxnBounds(Clostridium_symbiosum_WAL_14673,'EX_o2(e)',0,'l'); % ID: 1364
Bacteroides_thetaiotaomicron_VPI_5482=changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482,'EX_o2(e)',0,'l'); % ID: 535

models = {Clostridium_symbiosum_WAL_14673; Bacteroides_thetaiotaomicron_VPI_5482};
Names = {'Clostridium_symbiosum_WAL_14673'; 'Bacteroides_thetaiotaomicron_VPI_5482'};

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
crecimiento = zeros(length(models),2);

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

        if i == 1
            Clostridium_symbiosum_WAL_14673 = changeRxnBounds(Clostridium_symbiosum_WAL_14673, Clostridium_symbiosum_WAL_14673.rxns(cellfun(@isempty, strfind(Clostridium_symbiosum_WAL_14673.rxns, 'EX_'))==0), 0, 'l');
            Clostridium_symbiosum_WAL_14673 = changeRxnBounds(Clostridium_symbiosum_WAL_14673, medium, -rates, 'l');
            fba_model_Clostridium_symbiosum_WAL_14673= optimizeCbModel(Clostridium_symbiosum_WAL_14673, 'max');
            exc_idx_Clostridium_symbiosum_WAL_14673 = cellfun(@isempty, strfind(Clostridium_symbiosum_WAL_14673.rxns, 'EX_'))==0;
        end
        
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, Bacteroides_thetaiotaomicron_VPI_5482.rxns(cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0), 0, 'l');
        Bacteroides_thetaiotaomicron_VPI_5482 = changeRxnBounds(Bacteroides_thetaiotaomicron_VPI_5482, medium, -rates, 'l');
        
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482= optimizeCbModel(Bacteroides_thetaiotaomicron_VPI_5482, 'max');
        exc_idx_Bacteroides_thetaiotaomicron_VPI_5482 = cellfun(@isempty, strfind(Bacteroides_thetaiotaomicron_VPI_5482.rxns, 'EX_'))==0; 
        
        fba_model= optimizeCbModel(model_test);
        exc_idx = cellfun(@isempty, strfind(model_test.rxns, 'EX_'))==0;
                
        Clostridium_symbiosum_WAL_14673_mono_table = table(Clostridium_symbiosum_WAL_14673.rxns(find(fba_model_Clostridium_symbiosum_WAL_14673.v(:,1) ~= 0),1), ...
        fba_model_Clostridium_symbiosum_WAL_14673.v(find(fba_model_Clostridium_symbiosum_WAL_14673.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
      
        Bacteroides_thetaiotaomicron_VPI_5482_mono_table = table(Bacteroides_thetaiotaomicron_VPI_5482.rxns(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ...
        fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(find(fba_model_Bacteroides_thetaiotaomicron_VPI_5482.v(:,1) ~= 0),1), ... 
        'VariableNames',{'rxn','v'});
        
        writetable(Clostridium_symbiosum_WAL_14673_mono_table,'Clostridium_symbiosum_WAL_14673_monocultivo_agora.xlsx');
        writetable(Bacteroides_thetaiotaomicron_VPI_5482_mono_table,'Bacteroides_thetaiotaomicron_VPI_5482_monocultivo_agora.xlsx');

        crecimiento(i, j) = fba_model.f;
        sol_mono{end+1} = fba_model.x;
    end   
end

%% Graphical fluxes of external reactions/metabolites of interest in monocultures
targets = {'EX_inulin(e)'; 'EX_biomass(e)'; 'EX_ac(e)'; 'EX_but(e)'; 'EX_for(e)'; 'EX_lac_D(e)'; 'EX_lac_L(e)'; 'EX_succ(e)'; 'EX_ppa(e)'};
ids_Clostridium_symbiosum_WAL_14673 = findRxnIDs(Clostridium_symbiosum_WAL_14673, targets);
targets_ids_Clostridium_symbiosum_WAL_14673 = ids_Clostridium_symbiosum_WAL_14673(ids_Clostridium_symbiosum_WAL_14673 > 0);
targets_filt_Clostridium_symbiosum_WAL_14673 = targets(ids_Clostridium_symbiosum_WAL_14673 > 0);

ids_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(Bacteroides_thetaiotaomicron_VPI_5482, targets);
targets_ids_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Bacteroides_thetaiotaomicron_VPI_5482(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Bacteroides_thetaiotaomicron_VPI_5482 = targets(ids_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Clostridium_symbiosum_WAL_14673, '_', '\_')), fba_model_Clostridium_symbiosum_WAL_14673.x(targets_ids_Clostridium_symbiosum_WAL_14673));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')})
title('mZMB + Inulin Mono-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), fba_model_Bacteroides_thetaiotaomicron_VPI_5482.x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')})
title('mZMB + Inulin Mono-culture')
hold off

metabolitos_mono_Clostridium_symbiosum_WAL_14673_table = table(targets_filt_Clostridium_symbiosum_WAL_14673, ...
    fba_model_Clostridium_symbiosum_WAL_14673.x(targets_ids_Clostridium_symbiosum_WAL_14673), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_mono_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Bacteroides_thetaiotaomicron_VPI_5482, ...
    fba_model_Bacteroides_thetaiotaomicron_VPI_5482.x(targets_ids_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

%% Creation of the steadycom community model
EcCom = createMultipleSpeciesModel(models, Names);
% This is because the optimizer needs the csense field
EcCom.csense = char('E' * ones(1,numel(EcCom.mets))); 
% With this we will obtain the rates and reactions of each species, and those that are exclusive of the space between bacteria (lumen).
[EcCom.infoCom, EcCom.indCom] = getMultiSpeciesModelId(EcCom,Names);

%% Add biomass reaction information
rxnBiomass = EcCom.rxns(find(contains(EcCom.rxns, strcat(Names, 'biomass')))) % MODELOS AGORA
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

%% Optimization
options = struct();
options.GRguess = 0.1; % Initial value for growth rate
options.GRtol = 1e-4; % Tolerance for growth rate 
options.algorithm = 1; %  default
tic
[sol, result] = SteadyCom(EcCom, options);
toc

Clostridium_symbiosum_WAL_14673_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,1) ~= 0),1)), ... 
    'VariableNames',{'rxn','v'});

Bacteroides_thetaiotaomicron_VPI_5482_table = table(EcCom.infoCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2), ... 
    sol.full(EcCom.indCom.EXsp(find(EcCom.indCom.EXsp(:,2) ~= 0),2)), ... 
    'VariableNames',{'rxn','v'});

com_table = table(EcCom.rxns(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    sol.full(cellfun(@isempty, strfind(EcCom.rxns, 'EX_'))==0), ... 
    'VariableNames',{'rxn','v'});

%% export excel tables with consumed/produced metabolites with steadycomFBA
Clostridium_symbiosum_WAL_14673_table(~Clostridium_symbiosum_WAL_14673_table.v,:) = []; % https://www.mathworks.com/matlabcentral/answers/407401-delete-zero-rows-from-a-table
Bacteroides_thetaiotaomicron_VPI_5482_table(~Bacteroides_thetaiotaomicron_VPI_5482_table.v,:) = [];
com_table(~com_table.v,:) = [];

Clostridium_symbiosum_WAL_14673_table.rxn = regexprep(Clostridium_symbiosum_WAL_14673_table.rxn, '[u]tr', '(e)'); % me reemplza los _e por TT
Bacteroides_thetaiotaomicron_VPI_5482_table.rxn = regexprep(Bacteroides_thetaiotaomicron_VPI_5482_table.rxn, '[u]tr', '(e)'); % me reemplza los _e por TT
com_table.rxn = regexprep(com_table.rxn, '[u]tr', '(e)'); % me reemplza los _e por TT

writetable(Clostridium_symbiosum_WAL_14673_table,'Clostridium_symbiosum_WAL_14673_consorcio_con_Bacteroides_thetaiotaomicron_VPI_5482_agora.xlsx');
writetable(Bacteroides_thetaiotaomicron_VPI_5482_table,'Bacteroides_thetaiotaomicron_VPI_5482_consorcio_con_Clostridium_symbiosum_WAL_14673_agora.xlsx');
writetable(com_table,'EX_comunidad_Clostridium_symbiosum_WAL_14673_Bacteroides_thetaiotaomicron_VPI_5482_agora.xlsx');

%% Graph for community external metabolites such as butyrate
targets_com_Clostridium_symbiosum_WAL_14673 = {'Clostridium_symbiosum_WAL_14673IEX_inulin[u]tr';
    'Clostridium_symbiosum_WAL_14673IEX_kestopt[u]tr';
    'Clostridium_symbiosum_WAL_14673IEX_biomass[c]tr';
    'Clostridium_symbiosum_WAL_14673IEX_ac[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_but[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_for[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_lac_D[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_lac_L[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_succ[u]tr'; 
    'Clostridium_symbiosum_WAL_14673IEX_ppa[u]tr'}; 
    
 targets_com_Bacteroides_thetaiotaomicron_VPI_5482 = {'Bacteroides_thetaiotaomicron_VPI_5482IEX_inulin[u]tr';
     'Bacteroides_thetaiotaomicron_VPI_5482IEX_kestopt[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_biomass[c]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_ac[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_but[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_for[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_lac_D[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_lac_L[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_succ[u]tr'; 
    'Bacteroides_thetaiotaomicron_VPI_5482IEX_ppa[u]tr'};

targets_com = {'EX_biomass[c]';
    'EX_inulin[u]';
    'kestopt[e]';
    'EX_ac[u]'; 
    'EX_but[u]'; 
    'EX_for[u]'; 
    'EX_lac_D[u]'; 
    'EX_lac_L[u]'; 
    'EX_succ[u]'; 
    'EX_ppa[u]'};

ids_Com_Clostridium_symbiosum_WAL_14673 = findRxnIDs(EcCom, targets_com_Clostridium_symbiosum_WAL_14673);
targets_ids_Com_Clostridium_symbiosum_WAL_14673 = ids_Com_Clostridium_symbiosum_WAL_14673(ids_Com_Clostridium_symbiosum_WAL_14673 > 0);
targets_filt_Com_Clostridium_symbiosum_WAL_14673 = targets_com_Clostridium_symbiosum_WAL_14673(ids_Com_Clostridium_symbiosum_WAL_14673 > 0);

ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = findRxnIDs(EcCom, targets_com_Bacteroides_thetaiotaomicron_VPI_5482);
targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 = ids_Com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);
targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482 = targets_com_Bacteroides_thetaiotaomicron_VPI_5482(ids_Com_Bacteroides_thetaiotaomicron_VPI_5482 > 0);

ids_Com = findRxnIDs(EcCom, targets_com);
targets_ids_Com = ids_Com(ids_Com > 0);
targets_filt_Com = targets_com(ids_Com > 0);

figure
hold on
bar(categorical(strrep(targets_filt_Com_Clostridium_symbiosum_WAL_14673, '_', '\_')), sol.full(targets_ids_Com_Clostridium_symbiosum_WAL_14673));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, '_', '\_')), sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{2}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

figure
hold on
bar(categorical(strrep(targets_filt_Com, '_', '\_')), sol.full(targets_ids_Com));
ylabel('mmol gDW-1 h-1');
xlabel({strrep(Names{1}, '_', '\_') strrep(Names{2}, '_', '\_')});
title('mZMB + Inulin Co-culture')
hold off

metabolitos_com_Clostridium_symbiosum_WAL_14673_table = table(targets_filt_Com_Clostridium_symbiosum_WAL_14673, ...
    sol.full(targets_ids_Com_Clostridium_symbiosum_WAL_14673), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_com_Bacteroides_thetaiotaomicron_VPI_5482_table = table(targets_filt_Com_Bacteroides_thetaiotaomicron_VPI_5482, ...
    sol.full(targets_ids_Com_Bacteroides_thetaiotaomicron_VPI_5482), ... 
    'VariableNames',{'Metabolite','v'});

metabolitos_comunidad_totales_table = table(targets_filt_Com, ...
    sol.full(targets_ids_Com), ... 
    'VariableNames',{'Metabolite','v'});

%% SteadyCom-FVA
% maximum percentage of total biomass 
options.optBMpercent = 100;
options.rxnNameList = strcat('X_', Names); % rset the variables to look for relative abundances in the steadycom FVA
options.optGRpercent = [0:1:100]; % This option is to determine the percentages of the maximum growth rate that the FVA considers, thus how many FVA it makes, (in this case 0, 1, 2, ..., 100), you can change the range
%tolerance for final growth rate
options.GRtol = 1e-6;
options.threads = 0;

[fvaComMin, fvaComMax, minFD, maxFD, GRvector, result, LP] = SteadyComFVA(EcCom, options);

%% Envelope plot for relative abundance by community growth rate
grComV = result.GRmax * options.optGRpercent / 100; % Growth rates to be evaluated in the FVA
lgLabel = {strrep(Names{1}, '_', '\_');strrep(Names{2}, '_', '\_')};  % put all models

x_lower = 10;
x_upper = 0;
close all
figure
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(Names) % from 1 to number of mutant models 
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance'); % Relative abundance Bacteria growth rate (h^{-1})
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');
hold off

%% FVA optimizing for metabolite of interest (butyrate)
options.optBMpercent = 100;
options.rxnNameList = EcCom.rxns(find(contains(EcCom.rxns, 'IEX_but[u]')));
options.optGRpercent = [0:1:100];
options.GRtol = 1e-6;
[fvaComMin_B, fvaComMax_B, minFD_B, maxFD_B, GRvector_B, result_B, LP_B] = SteadyComFVA(EcCom, options);

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
    y_max=max(y);
    y = y(~isnan(y));
    if crecimiento(j,1) == max(crecimiento(:,1))
        if min(y) < x_lower 
            x_lower = min(y);
        end 
        if max(y) >= x_upper
            x_upper = y_max;
        end 
        if max(y) >= y_max
            y_max = max(y)
        end 
    end 
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

tl(1) = title('\underline{SteadyComFVA butyrate}', 'Interpreter', 'latex');
lg = legend(lgLabel);
% lg.Box = 'off';
yl(1) = ylabel('Butyrate flux (mmol gDW-1 h-1)');
xl(1) = xlabel('Community growth rate (h^{-1})');
hold off

%% Combined figures
lgLabel1 = {strrep(Names{1}, '_', '\_'); strrep(Names{2}, '_', '\_')}; % poner todos los modelos mutantes

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
for j = 1:length(Names)
    y = [fvaComMin(j, :), fliplr(fvaComMax(j, :))];
    y = y(~isnan(y));
    p(j, 1) = plot(x(~isnan(y)), y(~isnan(y)), 'LineWidth', 2);
end

legend(lgLabel1)
xlabel('Growth rate (h^{-1})');
ylabel('Relative abundance');
title('\underline{SteadyComFVA biomass}', 'Interpreter', 'latex');

% Plot the second subplot
subplot('Position', subplotPositions(2, :))
hold on
grid on
x = [grComV(:); flipud(grComV(:))];
y_max = 0;
for j = 1:length(options.rxnNameList)
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
lgLabel2 = {strrep(Names{1}, '_', '\_')}; % put butyrate producing bacteria here
lg = legend(lgLabel2);
yl(1) = ylabel('Butyrate flux (mmol gDW^{-1} h^{-1})');
xl(1) = xlabel('Community growth rate (h^{-1})');

% Save the combined figure
saveas(p(j, 1),'EX_combinado_Clostridium_symbiosum_WAL_14673_Bacteroides_thetaiotaomicron_VPI_5482_agora.png');
hold off

%% Monte Carlo Random Flux-sampling for the community model
% Required! Biomass restrictions before running sampling (steadycom-fba outputs were used):
EcCom = changeRxnBounds(EcCom,'Clostridium_symbiosum_WAL_14673biomass163',result.vBM(1),'b');
EcCom = changeRxnBounds(EcCom,'Clostridium_symbiosum_WAL_14673IEX_biomass[c]tr',result.vBM(1),'b');
EcCom = changeRxnBounds(EcCom,'Bacteroides_thetaiotaomicron_VPI_5482biomass345',result.vBM(2),'b');
EcCom = changeRxnBounds(EcCom,'Bacteroides_thetaiotaomicron_VPI_5482IEX_biomass[c]tr',result.vBM(2),'b');
EcCom = changeRxnBounds(EcCom,'EX_biomass[c]',result.GRmax,'b');

changeCobraSolverParams('LP', 'feasTol', 1e-4);

options.nStepsPerPoint = 100;
options.nPointsReturned=200000;
   
[modelSampling,samples]=sampleCbModel(EcCom, 'ModelSampling',  'CHRR', options);
 
%% Export sampling results
writematrix(samples,'samples_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_200000_vmh_ids.csv') 

%% Export sampling reaction names
EcCom_rxns = EcCom.rxns
writecell(EcCom_rxns,'samples_EcComRXNS_Bifidobacterium_animalis_lactis_PT33_Clostridium_innocuum_HFG2_vmh_ids.csv') 
