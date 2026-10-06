function get_nisBiosynthesisMac2
paths = soilSynComPaths();
% Load the final M9 TFA collection.
cd(paths.models.m9TFA)

% here you can give a list of GEM names
%names = {'Gilliamella_apicola','Snodgrassella_alvi'};

% or read them automatically from the folder
filenames = dir;
for l = 1:length(filenames)
    names{l,1} = filenames(l).name;
end
% make sure to keep only mat files
names = names(contains(names,'.mat'));
names = strrep(names,'.mat','');
cd(fileparts(mfilename('fullpath')))
addpath(genpath(fileparts(mfilename('fullpath'))))

for p = 15
 %check where it is saved
 paramNISbiosynthesis    = struct (...%model parameters
        'Organism',                           names{p}                           ,... %  strrep(model.description,'.xml','') human, ecoli, putida
        'GEMname',                            names{p}                              ,... % name of the GEM used (as it is saved in GEMs folder)
        'RedModelName',                        names{p}                 ,... % strrep(model.description,'.xml','') choose a name for the reduction
        'SelectedSubsystems',                 {{'Glycolysis / Gluconeogenesis';
                                                'Pentose phosphate pathway';
                                                'Citrate cycle (TCA cycle)';
                                                'Pyruvate Metabolism';
                                                'Glyoxylate and dicarboxylate metabolism'
                                                
                                               
                                                
                                                
                                              }}            ,... % default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'x';'y';'z'}}
        'AddETCAsSubsystem',                  'no'                                  ,... % yes, no
        'AddExtracellularSubsystem',          'no'                  ,... % no, default (default is defined in the organism&GEM case-file), customly defined in a cell e.g. {{'EX_x_e';'EX_y_e';'EX_z_e'}}
        'AerobicAnaerobic',                   'aerobic'                              ,... % aerobic, anaerobic
        'ListForInorganicMets',               'automatic'                            ,... % curated, automatic
        'ListForCofactorPairs',               'curated'                              ,... % curated, automatic
        'ZeroZeroGEMbounds',                  'Original'                             ,... % Original, DefineCustom, OpenTo100
        'case_filename',                      'case_soil2'                           ,... % name of the matlab function for the specific organism
        ...%redGEM parameters
        'L',                                  0                                      ,... % 1,2,...
        'D',                                  0                                      ,...
        'startFromMin',                       'yes'                                   ,... % yes, no
        'ThrowErrorOnDViolation',             'error'                                ,... % error/continue If any two subsystems cannont connect to the level D
        'OnlyConnectExclusiveMets',           'no'                                   ,... % yes, no
        'ConnectIntracellularSubsystems',     'yes'                                  ,... % yes, no
        'ApplyShortestDistanceOfSubsystems',  'bothways'                             ,... % bothways,eachdirection'
        ...%redGEMX parameters
        'performREDGEMX',                     'no'                                   ,... % yes, no
        'NumOfConnections',                   'SminMetE'                           ,... % OnePerMetE ,SminMetE
        ...%lumpGEM parameters
        'performLUMPGEM',                     'yes'                                  ,... % yes, no (do we want to perform lumping or not?)
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
        'CPLEX_PATH',                         paths.external.CPLEX                    ,... %provide path to CPLEX
        'TFA_PATH',                           paths.external.matTFA                   ,... %provide path to matTFA
        'thermo_data_PATH',                   paths.nisThermoData                     ,...% Provide the path to the thermodynamic data (including the file name)
        'redGEM_PATH',                        paths.external.redGEM                   ,... %provide path to redGEM
        'nis_PATH',                           paths.external.NIS                      ,... %provide path to NIS
        'output_PATH',                        fullfile(paths.outputRoot, 'NIS', [names{p} '_D0_Biosynthesis']));             % provide output path where the files will be saved
     nisBiosynthesis(paramNISbiosynthesis);
end
