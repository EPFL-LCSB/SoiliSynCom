%% Setting up the paths, solver and preparing the model for ReMIND
clear
clc

% Provide path to folder where your results will be saved
saving_directory = './ime_output';
if ~isdir(saving_directory)
    mkdir(saving_directory)
end

%add paths
%ReMIND_directory =  '${REMIND_ROOT}/';
ReMIND_directory =  '${REMIND_ROOT}/';
mattfa_directory = '/yourPATH';
cplex_directory =  '/yourPATH';
addpath(genpath(ReMIND_directory));
addpath(genpath(mattfa_directory));
addpath(genpath(cplex_directory));
changeCobraSolver('cplex_direct');

% Tag true to apply phenomapping with TFA (thermodynamic constraints will
% be taken into account)
tagThermo = 0;

% Prepare model for phenomapping
%load model and thermo data
%modelDescription = 'myModel';
%modelPath = '${SOILSYNCOM_PROJECT_ROOT}/bbet_PTYG.mat';
modelPath = '${SOILSYNCOM_DATA_ROOT}/cuvk_PTYG.mat';
% modeldescription = strrep(modelPath,'/yourPATH/matReMIND/models/','');
% modeldescription = strrep(modeldescription,'.mat','');
modeldescription = 'cuvk';
load(modelPath)
% model = ttmodel_rfa;
model.description = modeldescription;
idx_exchange = find(findExcRxns(model));
growthRate = 0.16;
f = find(contains(model.varNames,'F_EXC_'));
model.var_ub(f) = 100;
%% set envirnoment
%Define the envirnoment you want to consider for the ReMIND workflow
% block uptakes froom demand and exchange reactions
% allow only the uptakes you want to consider for the analysis
idx_exchange = find(findExcRxns(model));
model.lb(idx_exchange) = 0;

%for example
inorganics = {'EX_cl_e';'EX_ca2_e';'EX_cobalt2_e';'EX_mobd_e';...
    'EX_cu2_e';'EX_fe2_e';'EX_fe3_e';'EX_h2o_e';'EX_k_e';'EX_mg2_e';...
    'EX_mn2_e';'EX_na1_e';'EX_pi_e';'EX_so4_e';...
    'EX_zn2_e';'EX_o2_e';'EX_h_e';'EX_nh4_e'};
f = find(ismember(model.rxns,inorganics));
model.lb(f) = -25;


% vitamins and carbon source
vitamins = {'EX_thm_e';'EX_btn_e';...
    'EX_ribflv_e';'EX_pydx_e';'EX_4abz_e';'EX_lipoate_e';...
    'EX_nac_e';'EX_fol_e';'EX_cbl1_e';'EX_pnto__R_e'};
f = find(ismember(model.rxns,vitamins));
model.lb(f) = -1;

carbon = {'EX_glc__D_e'};
f = find(ismember(model.rxns,carbon));
model.lb(f) = -10;

% Bound secretions to get realistic DiMEs
% To bound all secretions
model.ub(idx_exchange) = 25;

% comment
%here I see 4hbz secretion alternates with thm uptake so I block the secretion to reduce the alternatives
% also fe3 and fe2 alternate

% block sinks
f = find(contains(model.rxns,'sink_'));
model.ub(f) = 0;
f = find(ismember(model.rxns,{'EX_h2co3_e';'EX_co_e';'EX_cellb_e'}));
model.ub(f) = 0;
% if growth rate if known
growthRate = 0.26;
%%
if tagThermo ==  0
    [model, checkList, tagReady] = initDiMEsModel(model,'',tagThermo);
else
    % this is the thermo database included in the mattfa toolbox
    thermo_data_directory = '/yourPATH/mattfa/thermoDatabases/thermo_data.mat';
    thermo_data = load(thermo_data_directory);
    thermo_data = thermo_data.DB_AlbertyUpdate;
    [model, checkList, tagReady] = initDiMEsModel(model,thermo_data,tagThermo);
end
%% Save prepared model and get carbon source name
carbonSource = strrep(carbon,'EXC_','');
carbonSource = strrep(carbonSource,'_e','');
save(strcat(saving_directory,modeldescription,'_4ReMIND.mat'), 'model');