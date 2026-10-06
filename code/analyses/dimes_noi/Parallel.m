% .. Author:
% Anush Chiappino-Pepe 2018
% Evangelia Vayena 2023
paths = soilSynComPaths();
addSoilSynComDependencies(paths);
initCobraToolbox();
changeCobraSolver('cplex_direct');
saving_directory = [fullfile(paths.outputRoot, 'DiMEs', 'FullAlternatives') filesep];
load(fullfile(paths.dimes, 'GrowthRates.mat'));
load(fullfile(paths.dimes, 'MetabolomicsPool.mat'));
if ~isfolder(saving_directory)
    mkdir(saving_directory)
end
tagThermo=0
Models{1}=load(fullfile(paths.models.richTFA, 'ThermobbetPTYGVitconstrained.mat')).ttmodel_bbet;
Models{5}=load(fullfile(paths.models.richTFA, 'ThermobuolPTYGVitconstrained.mat')).ttmodel_buol;
Models{3}=load(fullfile(paths.models.richTFA, 'ThermocaujPTYGVitconstrained.mat')).ttmodel_cauj;
Models{4}=load(fullfile(paths.models.richTFA, 'ThermocohmPTYGVitconstrained.mat')).ttmodel_cohm;
Models{2}=load(fullfile(paths.models.richTFA, 'ThermocpiPTYGVitconstrained.mat')).ttmodel_cpi;
Models{6}=load(fullfile(paths.models.richTFA, 'ThermocuvkPTYGVitconstrained.mat')).ttmodel_cuvk;
Models{7}=load(fullfile(paths.models.richTFA, 'ThermocxyPTYGVitconstrained.mat')).ttmodel_cxy;
Models{8}=load(fullfile(paths.models.richTFA, 'ThermodriPTYGVitconstrained.mat')).ttmodel_dri;
Models{9}=load(fullfile(paths.models.richTFA, 'ThermoflkPTYGVitconstrained.mat')).ttmodel_flk;
Models{10}=load(fullfile(paths.models.richTFA, 'ThermolguPTYGVitconstrained.mat')).ttmodel_lgu;
Models{11}=load(fullfile(paths.models.richTFA, 'ThermolurPTYGVitconstrained.mat')).ttmodel_lur;
Models{12}=load(fullfile(paths.models.richTFA, 'ThermomgkPTYGVitconstrained.mat')).ttmodel_mgk;
Models{13}=load(fullfile(paths.models.richTFA, 'ThermomipPTYGVitconstrained.mat')).ttmodel_mip;
Models{14}=load(fullfile(paths.models.richTFA, 'ThermomlotPTYGVitconstrained.mat')).ttmodel_mlot;
Models{15}=load(fullfile(paths.models.richTFA, 'ThermopmorPTYGVitconstrained.mat')).ttmodel_pmor;
Models{16}=load(fullfile(paths.models.richTFA, 'ThermopscfPTYGVitconstrained.mat')).ttmodel_pscf;
Models{17}=load(fullfile(paths.models.richTFA, 'ThermopzuPTYGVitconstrained.mat')).ttmodel_pzu;
Models{18}=load(fullfile(paths.models.richTFA, 'ThermoracdPTYGVitconstrained.mat')).ttmodel_racd;
Models{19}=load(fullfile(paths.models.richTFA, 'ThermorfaPTYGVitconstrained.mat')).ttmodel_rfa;
Models{20}=load(fullfile(paths.models.richTFA, 'ThermotroPTYGVitconstrained.mat')).ttmodel_tro;
Models{21}=load(fullfile(paths.models.richTFA, 'ThermovpePTYGVitconstrained.mat')).ttmodel_vpe;
% 
% Models = {'ttmodel_bbet', 'ttmodel_buol', 'ttmodel_cauj', 'ttmodel_cohm', ...
%          'ttmodel_cpi', 'ttmodel_cuvk', 'ttmodel_cxy', 'ttmodel_dri', 'ttmodel_flk', ...
%          'ttmodel_lgu', 'ttmodel_lur', 'ttmodel_mgk', 'ttmodel_mlot', ...
%          'ttmodel_pmor', 'ttmodel_pscf', 'ttmodel_pzu', 'ttmodel_racd', ...
%          'ttmodel_rfa', 'ttmodel_tro', 'ttmodel_vpe', 'ttmodel_mip'};
% Initialize parallel pool with 55 workers
parpool('local', 70);
maxNumCompThreads(2);
metsIdx = cell(1, length(Models));
immOutput = cell(1, length(Models));

parfor m=5
    model = Models{m};
    modeldescription = model.id;
    model.description = modeldescription;
    %     idx_exchange = find(findExcRxns(model));
    growthRate =GrowthRates(m);
    %     model.ub(idx_exchange) = 25;
    f=find(contains(model.varNames,'R_EXC'));
    model.var_ub(f)=0;
    f=find(contains(model.rxns,'EXC'));
    model.lb(f)=0;
    for h=1:length(EXCMetabo)
        f=find(ismember(model.rxns,EXCMetabo(h)));
        model.lb(f)=EXClbFBA(h);
    end
    for h=1:length(EXCMetabo)
        f=find(ismember(model.varNames,strcat('R_',EXCMetabo(h))));
        model.var_ub(f)=EXClb(h);
    end
    [model, checkList, tagReady] = initDiMEsModel(model,'',tagThermo);
%     save(strcat(saving_directory,'DiMEsReadyModels/',modeldescription,'_4ReMIND.mat'), 'model');
    
% FIRST RUN file_setup.m for each one of the models and then run this script for all the models
    %%%%%%%%%%%%%%%%%%%%%%
    % DiMEs ANALYSIS %
    %%%%%%%%%%%%%%%%%%%%%%
    % inputs
    carbonSource = {''};
    sol=solveTFAmodelCplex(model);
    essThr = 0.99; % essentiality threshold (% of optimal growth to be defined, see minObj). If a knockout leads to grow below this threshold the gene will be considered essential for growth.
    minObj = essThr*growthRate; %essThr*grRate; % minimal required growth
    maxObj = growthRate; % upper bound in growth
    drainsForiMM = {}; % names of drains/exchanges used for IMM analysis. Empty means by default it will search for minimal media accross all drains/exchanges
    metabData = []; % provide metabolomics data to integrate in this analysis. Empty if none.
    NumAlt = 3000; % number of alternative minimal media to find. We define 1 for you to test that the pipeline works. But it is suggested to define 5000 or more to make sure you identify all alternatives
    time = 600; % time limit for optimization in seconds. If empty we do not set a time limit.
    tagMin = 0; % additional constrain to avoid generating suboptimal solutions, meaning that if 1 we will not identify media that is not minimal
    tagBoth = 1; % 1 --> iMEs , 0 --> iMMs
    tagOptYield = 1; % adding parsimonious uptake constraint
    weight = 'CA'; %weights on parsimonious constraint 'MW' (molecular weight),'CA' (carbon atoms),'default'(all w=1)
    tagCUT_pLOW=1; %add a lower bound on uptakes, this is used to constrain the solution space within a yiled cut
    tagMinSize=1; %add the block condition so that size of the alternatives <= min_size+2
    filename = strcat(modeldescription,'_PhenoMappingExchanges');
    p= [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
    p= 1./p; % percentages for parsimonious, here you can discretise as you wish
    
    
    EXCMetaboIdx=find(ismember(model.rxns,EXCMetabo));
    EXCMetaboIdxTFA=find(ismember(model.varNames,strcat('R_',EXCMetabo)));
    % all carbon sources
    [idx,carbonExchange] = getCarbonSources(model,'tfa'); %gets all carbon uptake reactions
    carbonExchange=intersect(carbonExchange,EXCMetaboIdx);
    idx=intersect(idx,EXCMetaboIdxTFA);
    metsIdx{m} = zeros(length(carbonExchange), 1);
    for i = 1:length(carbonExchange)
        metsIdx{m}(i) = find(model.S(:, carbonExchange(i)));
    end
    
    %Weights on parsiminious constraint
    if strcmp(weight,'CA') % get carbon atoms
        w = getCarbonAtoms(model,metsIdx{m});
    elseif strcmp(weight,'MW') % get MW(molecular weights)
        w = computeMW(model,model.mets(metsIdx{m}));
    elseif strcmp(weight,'default') %set weights as
        w = ones(length(idx),1);
    end
    
    %Constraint growth (f=objectve function)
    model.var_ub(find(model.f)) = maxObj;
    model.var_lb(find(model.f)) = minObj;
    
    idxTFA = getCarbonSources(model,'tfa');
    idxTFA=intersect(idxTFA,EXCMetaboIdxTFA);
    idx_C = find(ismember(model.rxns,strrep(model.varNames(idxTFA),'R_','')));
    
    % Parsimonious uptake problem
    model.f = zeros(length(model.varNames),1);
    model.f(idx) = w; %set the weights for the exchange reactions
    
    %Minimize
    model.objtype = 1; % 1-> minmize , -1 -> maximize
    sol_pUP = solveTFAmodelCplex(model,time); %solve to minimize flux of exchange reactions
    
    model.f = zeros(length(model.varNames),1);
    %Maximize
    model.objtype = -1;
    model.f(find(ismember(model.varNames,strcat('F_',model.rxns(find(model.c)))))) = 1; %put objective back to maximize growth
    
    
    %not binary variables for the inorganics
    
    
    [model_p,drains,min_size] = analysisDiMEs(model,tagBoth,minObj, maxObj, drainsForiMM, ...
        metabData,time);
    changeCobraSolver('cplex_direct');
    % Parsimonious uptake constraint
    [NumCons, NumVars] = size(model_p.A);
    NewCons = zeros(NumVars,1);
    NewCons(idx) = w;
    model_p.A(NumCons+1,:) = NewCons;
    model_p.rhs(NumCons+1) = 1.01*p(1)*sol_pUP.val;
    model_p.constraintNames{NumCons+1} = ['CUT_pUP'];
    model_p.constraintType{NumCons+1} = '<';
    
    if tagCUT_pLOW
        [NumCons, NumVars] = size(model_p.A);
        NewCons = zeros(NumVars,1);
        NewCons(idx) = w;
        model_p.A(NumCons+1,:) = NewCons;
        model_p.rhs(NumCons+1) = 0.99*p(1)*sol_pUP.val;
        model_p.constraintNames{NumCons+1} = ['CUT_pLOW'];
        model_p.constraintType{NumCons+1} = '>';
    end
    
    %here we can run sequentially for the different yield cuts but also in
    %parallel (i.e. different matlab instances)
    for i = 10:length(p)
        f = find(ismember(model_p.constraintNames,{'CUT_pUP'}));
        model_p.rhs(f) = 1.01*p(i)*sol_pUP.val;
        f = find(ismember(model_p.constraintNames,{'CUT_pLOW'}));
        model_p.rhs(f) = 0.99*p(i)*sol_pUP.val;
        
        
        %printConstraint(model_p,{'CUT_pUP'})
        fprintf('Generating alternatives for %4.2f % of the optimal biomass yield \n',1/p(i)*100)
        saving_directory_new = strcat(strcat(saving_directory,num2str(1/p(i)*100),'%_opt_yield_',cell2str(carbonSource),'/'));
        modeldescription=model_p.id;
        mkdir(saving_directory_new);
        
        sol = solveTFAmodelCplex(model_p,time);
        if ~isempty(sol.val)
            
            [DPsimm,modelmilp2] = findDPMax_DiMEs(model_p, NumAlt, model_p.indUSE,time, tagMin, strcat(saving_directory_new,filename),tagMinSize);
            r1=saver1(saving_directory_new,filename,model_p,DPsimm);
%             % save workspace
%             save(strcat(saving_directory_new,filename,'.mat'));
%             
%             % save results
%             filename = strcat(modeldescription,'_PhenoMappingExchanges');
%             r1 = load(strcat(saving_directory_new,filename,'.mat'));
            
            if not(isempty(DPsimm))
                % Extract info about composition of the IMMs
                [immOutput{m}.StatsMets, immOutput{m}.Mets, immOutput{m}.drainClass, ...
                    immOutput{m}.statsClass,immOutput{m}.result,immOutput{m}.Lengths,immOutput{m}.Alt,Active] = extractInfo_DiMEs(r1.model_p, r1.DPsimm, ...
                    r1.model_p.indUSE,  r1.model_p.rowDrain_F,...
                    r1.model_p.rowDrain_R,r1.model_p.F_Flux,  r1.model_p.R_Flux) ;
                
                
                immOutput{m}.summary = [strrep(immOutput{m}.Mets(:,1),',',''), ...
                    num2cell(immOutput{m}.result),num2cell(immOutput{m}.Lengths),num2cell(immOutput{m}.Alt),num2cell(immOutput{m}.StatsMets)];
                
                writeData(strcat(saving_directory_new,filename,'_',num2str(1/p(i)*100),'.csv'), immOutput{m}.summary,...
                    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
                    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
                    '%s\t%s\t%s\t%s\t%s');
%                 saveDiMEsData(saving_directory_new, modeldescription)
%                 save(strcat(saving_directory_new,'Active_',modeldescription,'.mat'),'Active');
%                 clear Active
                saveDiMEs(saving_directory_new, model_p,Active)
                elimIntermPMFile(saving_directory_new,filename)
            end
            f = find(contains(model_p.constraintNames,'CUTimm'));
            model_p.constraintNames(f) = strcat(model_p.constraintNames(f),'_old');
        end
    end
    modeldescription=model_p.id;
    % Get unique solutions and prepare the input for the ILP step of ReMIND
    getUniqueSol(p,modeldescription,saving_directory)
    generateSummaryTables(model_p)
    delete(gcp('nocreate'))
end

function r1=saver1(saving_directory_new,filename,model_p,DPsimm)
 % save workspace
            modeldescription=model_p.id;
            save(strcat(saving_directory_new,filename,'.mat'),'model_p','DPsimm');
            % save results
            filename = strcat(modeldescription,'_PhenoMappingExchanges');
            r1 = load(strcat(saving_directory_new,filename,'.mat'));
end
function saveDiMEs(saving_directory_new, model_p,Active)
                modeldescription=model_p.id;
                save(strcat(saving_directory_new,'Active_',modeldescription,'.mat'),'Active');
end
function generateSummaryTables(model_p)
%Define the parsimonious values we used
modeldescription=model_p.id;
paths = soilSynComPaths();

rowsP=[];
rowsS=[];
rownames={};
saving_directory = [fullfile(paths.outputRoot, 'DiMEs', 'FullAlternatives') filesep];
load(strcat(saving_directory,modeldescription,'_unique_active.mat'))
Drains = Drains';
for i = 1:size(Active,2)
    max_P(i) = length(find(Active(:,i)==1));
end
for i = 1:size(Active,2)
    max_S(i) = length(find(Active(:,i)==-1));
end
max_P = max(max_P);
max_S= max(max_S);
P = cell(max_P,1)';
S = cell(max_S,1)';

for j=1:size(Active,2)
    idx_p_j=find(Active(:,j)==1);
    idx_s_j=find(Active(:,j)==-1);
    newrowsP=[Drains(idx_p_j)];
    newrowsS=[Drains(idx_s_j)];
    P(j,1:length(idx_p_j)) = newrowsP;
    S (j,1:length(idx_s_j)) = newrowsS;
    rownames= horzcat(rownames,strcat(modeldescription,'_l_',num2str(tag(j)),'_alt_',num2str(j)));
    newrowsP={};
    newrowsS={};
end


columnNamesP={};
columnNamesS={};
for i= 1:size(P,2)
    columnNamesP{i}=strcat('P_',int2str(i));
end
for i= 1:size(S,2)
    columnNamesS{i}=strcat('S_',int2str(i));
end
Allrows=horzcat(P,S);
idx = find(cellfun(@isempty, Allrows));
Allrows(idx)={'none'};
Allrows=regexprep(Allrows,"'","");
AllColumnNames=horzcat(columnNamesP,columnNamesS);
SummaryTable=array2table(Allrows);
SummaryTable.Properties.VariableNames=AllColumnNames;
SummaryTable.Properties.RowNames=rownames;
save(strcat(saving_directory,modeldescription,'_summary.mat'),'SummaryTable')
writetable(SummaryTable,strcat(saving_directory,modeldescription,'_summary.csv'),'WriteRowNames',true);
end
