function get_nisExtracellular
% go to nis/runFiles to run this function
cd ../
cd ./GEMs

% here you can give a list of GEM names
%names = {'Snodgrassella_alvi','Gilliamella_apicola'};

% or read them automatically from the folder
filenames = dir;
for i = 1:length(filenames)
    names{i,1} = filenames(i).name;
end
% make sure to keep only mat files
names = names(contains(names,'.mat'));
names = strrep(names,'.mat','');
cd ../
addpath(genpath('./'))

for i = 1: length(names)
 %check where it is saved
 paramNISextra    = struct (...%model parameters
        'Organism',                           names{i}                           ,... %  strrep(model.description,'.xml','') human, ecoli, putida
        'GEMname',                            names{i}                              ,... % name of the GEM used (as it is saved in GEMs folder)
        'RedModelName',                        names{i}                 ,... % strrep(model.description,'.xml','') choose a name for the reduction
        'SelectedSubsystems',                 {{'Glycolysis/Gluconeogenesis';
                                                'Pentose Phosphate Pathway';
                                                'Citric Acid Cycle';
                                                'Pyruvate Metabolism'
                                              }}            ,... % default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'x';'y';'z'}}
        'AddETCAsSubsystem',                  'yes'                                  ,... % yes, no
        'AddExtracellularSubsystem',          'no'                  ,... % no, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'EX_x_e';'EX_y_e';'EX_z_e'}}
        'AerobicAnaerobic',                   'aerobic'                              ,... % aerobic, anaerobic
        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100
        'case_filename',                      'case_nis'                           ,... % name of the matlab function for the specific organism
        ...%redGEM parameters
        'L',                                  0                                      ,... % 1,2,...
        'D',                                  0                                      ,...
        'startFromMin',                       'yes'                                   ,... % yes, no
        'ThrowErrorOnDViolation',             'error'                                ,... % error/continue If any two subsystems cannont connect to the level D
        'OnlyConnectExclusiveMets',           'no'                                   ,... % yes, no
        'ConnectIntracellularSubsystems',     'yes'                                  ,... % yes, no
        'ApplyShortestDistanceOfSubsystems',  'bothways'                             ,... % bothways,eachdirection'
        ...%redGEMX parameters
        'performREDGEMX',                     'yes'                                   ,... % yes, no
        'NumOfConnections',                   'SminMetE'                           ,... % OnePerMetE ,SminMetE
        ...%lumpGEM parameters
        'performLUMPGEM',                     'no'                                  ,... % yes, no (do we want to perform lumping or not?)
        'PercentOfmuMaxForLumping',            100                                   ,... % Please specify the percentage of muMax that we should impose for the lumped reactions (100, 90, ..., 10, 0)
        'addGAM',                             'no'                                  ,... % Would you like to extract and add to the reduced model a growth associated maintenance (GAM) reaction? yes, no
        'PreventBBBuptake',                   'no'                                   ,... % yes/no: if yes, do not allow uptake through any bbb drains.
        'NumOfLumped',                        'Smin'                                 ,... % OnePerBBB, Smin, Sminp1, Sminp2, Sminp3
        'ImposeThermodynamics',               'yes'                                  ,... % would you like to impose thermodynamic cnstraints? yes, no
        'AlignTransportsUsingMatFile',        'yesautomatic'                         ,... % yesusingmatfile, yesusingReducedmatfile, yesautomatic,no
        ...%postprocessing parameters
        'performPostProcessing',              'no'                                  ,... % yes, no
        ...%solver parameters
        'TimeLimitForSolver',                 'yes'                                  ,... % yes, no
        'CplexParameters',                    'LCSBDefault'                          ,... % CplexDefault, LCSBDefault
        ....%paths to folders
        'CPLEX_PATH',                         '/Applications/CPLEX_Studio1210/cplex/matlab'                        ,... %provide path to CPLEX
        'TFA_PATH',                           '${MATTFA_ROOT}'                      ,... %provide path to matTFA
        'thermo_data_PATH',                   strcat('./matFiles/thermo_data.mat')  ,...% Provide the path to the thermodynamic data (including the file name)
        'redGEM_PATH',                         '${REDGEM_ROOT}'                      ,... %provide path to redGEM
        'nis_PATH',                         '${NIS_ROOT}/'                      ,... %provide path to NIS
        'output_PATH',                        strcat('./output/',names{i},'_D0_EX'));                              % provide output path where the files will be saved
     nisExtracellular(paramNISextra);
end

