function nisBiosynthesis(varargin)
% This function generates subnetworks to produce biomass building blocks
% from a D = 0 core model of a GEM (GEnomescale Model).
% Adapted from the redGEM.m function
%
% INPUTS
% - L: is an adjustable parameter indicating how far from each subsystem is
%   the algorithm allowed to search in order to find connections among
%   the subsystems. More precicely, it indicates how distant, in terms of
%   reactions, can the metabolites of any two core sybsystems be from any
%   potential interconnecting metabolite.
% - D: given that the algorithm has searched paths of up to length L, D
%   specifies how many lengths to include in the new model. if not
%   specified, this is equal to L.
% - startFromMin: default = 'yes'. if 'yes', indicates that the L's that are
%   included in the model start at the shortest distance between subsystem
%   pairs. if zero, the model is built with the first D=L levels regardless
%   of the presence of paths.
% - viewStats: options to plot statistics on the connecting pathways.
%   default is 'no'
% - OnlyConnectExclusiveMets: if this is 'yes', than only metabolites that
%   are exlusive to each subsystem will be connected. if 'no', even if a
%   metabolite is shared, it will still be connected. default = 'yes'
% - ConnectIntracellularSubsystems: if this is 'yes', than metabolites that are
%   within a subsystem will be connected to eachother. if 'no', they will
%   not be connected. default is 'no';
%   ----------------------------------------------------------------------
%   Example: L=10 D=4 startFromMin= 'yes'
%   the adjacency matrix will be calculated up to length 10, but only the
%   first 4 lengths that connect subsystems will be included. If subsystem
%   A and B have a D=1 connection, then the paths of length 1:4 will be
%   included in the model. However, the shortest distance between
%   subsystems A and C is 3 so all paths of length 3:7 will be included.
%   ----------------------------------------------------------------------
% - stats: options to plot statistics on the connecting pathways. default=0

% OUTPUTS
% - RedModel: Reduced model
% - rxns_total:
% - activeRxns:
% - overall_am:
% - overall_am_names_second:
% - DPs_all:

tic

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SETTING UP THE OPTIONS OF THE REDUCTION  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set all non-specified flags to empty. In this way, whatever has not been
% specified by the "varargin", will be later on asked with individual
% questions. All the variables that have been specified in "varargin" will
% overwrite these empty settings.
[paramNISbiosynthesis, paramNames] = nisQuestionnaire(varargin{1});


% Assign the value names to all the parameters of the current redGEM run
eval(cell2mat(cellfun(@(x) [x,'= getfield(paramNISbiosynthesis,''',x,''');'],paramNames,'UniformOutput', false)'));

if ~exist([output_PATH,'/WorkSpaces/',Organism,'/',GEMname],'dir')
    mkdir([output_PATH,'/WorkSpaces/',Organism,'/',GEMname])
end
connectingpaths_folder = [output_PATH,'/ConnectingPaths'];
if ~exist(connectingpaths_folder,'dir')
    mkdir(connectingpaths_folder)
end
if ~exist([output_PATH,'/Models/',Organism],'dir')
    mkdir([output_PATH,'/Models/',Organism])
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SET NECESSARY THE PATHS AND CHECK THE VERSIONS USED  %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Restore the default path
restoredefaultpath

% Add required paths:
% - Cobra and the tFBA path
addpath(genpath(TFA_PATH))

% - Cplex path
addpath(genpath(CPLEX_PATH))

if isempty(which('cplex.p'))
    error('You need to add CPLEX to the Matlab-path!!')
end
% - NIS path
[path, ~, ~] = fileparts(mfilename('fullpath'));
addpath(genpath(path))

% - redGEM path
addpath(genpath(redGEM_PATH))

% - NIS path
addpath(genpath(nis_PATH))

% Set the solver
changeCobraSolver('cplex_direct','LP')
% Set properly the desired parameters for cplex LP and MILP
[mipTolInt, scalPar, feasTol, emphPar] = setCplexParamHelper(CplexParameters);

%% %%%%%%%%%%%%%%%%%%%%%%
%% LOAD NECESSARY FILES %
%% %%%%%%%%%%%%%%%%%%%%%%
% Load the thermodynamic database used
load(thermo_data_PATH)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set the properties of the GEM that is reduced %
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eval(['[OriginalGEM, GEMmodel, core_ss, Biomass_rxns, met_pairs_to_remove, InorgMetSEEDIDs, BBBsToExclude, ExtraCellSubsystem, OxPhosSubsystem,InorganicMetID] = '...
    case_filename,'(GEMname, ZeroZeroGEMbounds, ListForInorganicMets, ListForCofactorPairs, SelectedSubsystems, AddExtracellularSubsystem, DB_AlbertyUpdate);'])
%%
% Model Reaction Redundancy Check
redundant_table = ReacRedundancy(GEMmodel);
%% Prevent (close) the uptakes of biomass building blocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The reason why we have to do this is because otherwise the reduction
% might skip the generation of lumped reactions for some biomass building
% blocks.
GEMmodel_orig=GEMmodel;
n=find(ismember(GEMmodel.rxns,'ATPM'));
GEMmodel.lb(n)=0;
GEMmodel.ub(n)=0;
n=find(ismember(GEMmodel.varNames,'F_ATPM'));
GEMmodel.var_ub(n)=0;
n=find(ismember(GEMmodel.varNames,'R_ATPM'));
GEMmodel.var_ub(n)=0;
fprintf('disabled ATPM\n')

if strcmp(PreventBBBuptake, 'yes')
    GEMmodel = preventBBBuptake(GEMmodel);
end
%% Adding ETC as a subsystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(AddETCAsSubsystem,'yes')
    fprintf('The ETC has been added to the core subsystems\n')
    GEMmodel = generate_ETC_SubSyst(GEMmodel,OxPhosSubsystem);
    core_ss = [core_ss; {'ETC_Rxns'}];
elseif strcmp(AddETCAsSubsystem,'no')
    fprintf('The ETC has not been added to the core subsystems\n')
else
    error('Wrong option!')
end
%% Adding extracellular metabolites as a subsystem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(ExtraCellSubsystem)
    [ba] = contains(GEMmodel.rxns,ExtraCellSubsystem);
    GEMmodel.subSystems(ba) = {'ExtraCell'};
    core_ss = [core_ss; {'ExtraCell'}];
    fprintf('Generic ExtraCellular Metabolites have been added to the core subsystems\n')
elseif strcmp(AddExtracellularSubsystem,'no')
    fprintf('Did not add The ExtraCellular Metabolites as a core subsystem\n')
else
    error('Wrong option!')
end
%% For the process of lumping, we can use the entire genomescale model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GSM_ForLumping = GEMmodel;

if strcmp(ImposeThermodynamics, 'no')
    % Remove unnecessary thermo fields (if they exist). If not, functions
    % like the extractSubnetwork might have an error.
    GSM_ForLumping = removeThermoFields(GSM_ForLumping);
end
%% Connecting Subsystems
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If we only want to perform the lumping on a predefined subsystem, or set
% of subsystems (i.e. we do not want to connect these subsystem, and/or
% include more reactions) then we can directly jump to the lumping part.
% The above implies that D=0, as we do not wish to expand the

% To find the adjacency matrices and perform all the graph search to obtain
% the desired reaction network we need a genomescale model without biomass,
% and without the metabolite pairs:
GSM_ForAdjMat = GEMmodel;
% Find and remove the biomass reactions
[BMIds] = contains(GSM_ForAdjMat.rxns,Biomass_rxns);

GSM_ForAdjMat   = removeRxns(GSM_ForAdjMat, GSM_ForAdjMat.rxns(BMIds));

% We directly remove cofactors and small ionrganic metabolites from the S matrix
Smatrix = removeMetPairsFromS(GSM_ForAdjMat, met_pairs_to_remove, InorganicMetID);

GSM_ForAdjMat.S = Smatrix;
%% No Connections

rxns_ss = cell(length(core_ss),1);
mets_ss = cell(length(core_ss),1);
% For each of the subsystems:
%Omar Keshk
% Loop through each entry in model.subSystems
for i=1:length(core_ss)
    core = false(size(GSM_ForAdjMat.subSystems));
    for j = 1:length(GSM_ForAdjMat.subSystems)
        currentSubsystems = GSM_ForAdjMat.subSystems{j};
        
        % Loop through each element within the cell array if it's a cell array
        if iscell(currentSubsystems)
            for f = 1:length(currentSubsystems)
                % Compare each element to OxPhosSubsystem
                if strcmp(currentSubsystems{f}, core_ss{i})
                    core(j) = 1;
                    break; % If a match is found, no need to check further
                end
            end
        else
            % If it's not a cell array, compare directly
            core(j)=strcmp(currentSubsystems, core_ss{i});
        end
    end
    
    % - find the reactions of this subsystem, and their indices
    %         [~, rxns_ss{i}]=ismember(GSM_ForAdjMat.subSystems, core_ss{i});
    rxns_ss{i}=core;
    rxns_ss{i}=find(rxns_ss{i});
    mets_ss{i}=find(any(GSM_ForAdjMat.S(:,rxns_ss{i}),2));
    GSM_ForAdjMat.mets(mets_ss{i});
end

all_core_rxns = unique([vertcat(rxns_ss{:})]);

submodel = extractSubNetwork(GEMmodel, GSM_ForAdjMat.rxns(all_core_rxns));
% We go through all GS-reactions, and check if there exists a reactions
% that involves ONLY core metabolites. If yes, we keep it, and add it as a
% core (i.e. tightening)
k = 1;
for i=1:length(GEMmodel.rxns)
    mets = GEMmodel.mets(find(GEMmodel.S(:, i)));
    [~, ba] = ismember(mets, submodel.mets);
    if isempty(find(ba==0))
        core_rxns(k,1) = i;
        k=k+1;
    end
end
clear k;

otherReactions = setdiff(1:length(GSM_ForLumping.rxns),core_rxns);

% Transform indices to GEM model (GSM_ForLumping) indexing
[~, otherReactionsGSMForLump_idx] = ismember(GSM_ForLumping.rxns(otherReactions), GSM_ForLumping.rxns);

%% Prevent thermodynamic relaxation of ETC reactions and/or ATP-synthase
ATPsynth_RxnNames = [];
ATPsynth_RxnIds = [];
% In the case of thermodynamic feasibility when we add the thermodynaic
% constraints, there might exist reactions in the model for which we would
% prefer not to allow relaxation of their corresopnding DGo bounds. Such
% reactions are
% - ATP synthases, e.g.:
% 'adp_c + pi_c + 4 h_p <=> atp_c + 3 h_c + h2o_c'
% adp    'cpd00008'
% pi     'cpd00009'
% h      'cpd00067'
% atp    'cpd00002'
% h2o    'cpd00001'
if find(contains(GSM_ForLumping.metSEEDID,'cpd00008'))
SEEDIDs_atpSynth = {'cpd00008','cpd00009','cpd00067','cpd00002','cpd00001'}';
elseif find(contains(GSM_ForLumping.metSEEDID,'C00008'))
SEEDIDs_atpSynth = {'C00008','C00009','C00080','C00002','C00001'}';
end
% Stoich matrix with binary entries
Stoic01 = GSM_ForLumping.S~=0;
for i = 1:size(Stoic01, 2)
    SEEDIDs_Rxn_i = unique(GSM_ForLumping.metSEEDID(Stoic01(:,i)));
    if isempty(setxor(SEEDIDs_Rxn_i, SEEDIDs_atpSynth))
        ATPsynth_RxnIds = [ATPsynth_RxnIds i];
    end
end
ATPsynth_RxnNames = GSM_ForLumping.rxns(ATPsynth_RxnIds);
% - and reactions of the ETC sybsystem
%OmarKeshk
% ETC_RxnNames = GSM_ForLumping.rxns(find_cell('ETC_Rxns', GSM_ForLumping.subSystems));
ETC_RxnNames=[];
RxnNames_PrevThermRelax =  ATPsynth_RxnNames
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_Biosynthesis_0.mat;'])    %

[activeRxns, LumpedRxnFormulas, bbbNames, DPsAll, IdNCNTNER, relaxedDGoVarsValues_ForEveryLumpedRxn] = ...
        addBIOMASSnis(GSM_ForLumping, otherReactionsGSMForLump_idx, DB_AlbertyUpdate, BBBsToExclude, AerobicAnaerobic, ...
                   Organism, AlignTransportsUsingMatFile, TimeLimitForSolver, NumOfLumped, CplexParameters, ...
                   GEMname, RxnNames_PrevThermRelax, Biomass_rxns, ATPsynth_RxnNames, addGAM, PercentOfmuMaxForLumping, ...
                   ImposeThermodynamics, output_PATH);
    
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_Biosynthesis_a.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
    [~, indCORE] = ismember(1:length(GSM_ForLumping.rxns), IdNCNTNER);
    indCORE = find(indCORE==0);
    [RedModel, rxns ,info_LMDP] = ...
        add_lumped_reactions_nis(GSM_ForLumping, LumpedRxnFormulas, bbbNames, DB_AlbertyUpdate, ...
                             GSM_ForLumping.rxns(indCORE), ImposeThermodynamics, CplexParameters, Organism, GEMname, activeRxns, output_PATH);
   
    % > > > > > > > > > SAVING  WORKSPACE > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > > >
    [dateStr, timeStr] = getDateTimeStrings(date,clock);                                                      %
    eval(['save ',output_PATH,'/WorkSpaces/',Organism,'/',GEMname,'/',dateStr,'_',timeStr,'_',mfilename,'_Biosynthesis_b.mat;'])    %
    % < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < < <
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% ->->->->->->  PREPARE OUTPUT ->->->->-> %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Add date to the generation of the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[dateStr, timeStr] = getDateTimeStrings(date,clock);

% Create a model name based on the date and time of the generation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model_name = [RedModelName,'_',dateStr,'_',timeStr];
% Store the core reduced model stoichiometry in the output (without the lumped reactions and without the transports)
RedModel.core_model = extractSubNetwork(GEMmodel, GEMmodel.rxns(core_rxns));
%% Appending necessary informative fields to the model
RedModel.OriginalGEM = OriginalGEM;
RedModel.all_core_rxns = all_core_rxns;
RedModel.core_rxns = core_rxns;
% if L == 0
%     RedModel.mets_ss_ID = mets_ss_ID;
%     RedModel.rxns_ss_ID = rxns_ss_ID;
% end
% Store in the model some of the properties of the reduction
RedModel.ReductionProperties = paramNISbiosynthesis;

% Store in the model the core/initial subsystems that are selected for the
% reduction
RedModel.initialSubSystems = core_ss;

% Add a description field. ATTENTION: this is very important for the ORACLE
% workflow later on!!
RedModel.description = [RedModelName,'_Date_',dateStr,'_',timeStr,'_Generated_by_redGEMv1'];

% workflow later on!!
RedModel.solverPath = CPLEX_PATH;

% Add a field with the name of the user that generated this model
%user_str = regexp(pwd,'[\\/]Users[\\/](\w+)[\\/]','tokens');
user_str = 'Omar';
RedModel.GeneratedByUser = user_str;

%% Saving
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% How much time it took for all these?
rG_RT = toc;

% Save the runtime of the generation
RedModel.ReductionRuntime = rG_RT;

% Save the model in the output folder
eval([model_name,' = RedModel;']);
eval(['save ',output_PATH,'/Models/',Organism,'/',model_name,'.mat ', model_name,';']);
fprintf('nisBiosynthesis analysis completed')
end
