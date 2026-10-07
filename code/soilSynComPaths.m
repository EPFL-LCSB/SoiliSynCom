function paths = soilSynComPaths()
%SOILSYNCOMPATHS Return repository-relative data, model, and output paths.
%
% Add the repository code folder to the MATLAB path, then call:
%   paths = soilSynComPaths();
%
% External software locations are read from environment variables so the
% same scripts can be used on different workstations.

codeRoot = fileparts(mfilename('fullpath'));
repositoryRoot = fileparts(codeRoot);

paths.repositoryRoot = repositoryRoot;
paths.codeRoot = codeRoot;
paths.dataRoot = fullfile(repositoryRoot, 'data');
paths.workflowInputs = fullfile(paths.dataRoot, 'workflow_inputs');
paths.archiveRoot = getenv('SOILSYNCOM_ARCHIVE_ROOT');
paths.archiveWorkflowInputs = '';
if ~isempty(paths.archiveRoot)
    paths.archiveWorkflowInputs = fullfile(paths.archiveRoot, 'workflow_inputs');
end
paths.outputRoot = getenv('SOILSYNCOM_OUTPUT_ROOT');
if isempty(paths.outputRoot)
    paths.outputRoot = fullfile(repositoryRoot, 'outputs');
end

paths.thermodynamics = fullfile(paths.workflowInputs, 'thermodynamics');
paths.thermoData = fullfile(paths.thermodynamics, 'thermo_data.mat');
if ~isfile(paths.thermoData) && ~isempty(paths.archiveRoot)
    paths.thermoData = fullfile(paths.archiveWorkflowInputs, 'thermodynamics', 'thermo_data.mat');
end
paths.metaThermoData = fullfile(paths.thermodynamics, 'meta_thermo_data.mat');
paths.nisThermoData = fullfile(paths.thermodynamics, 'DB_AlbertyUpdate.mat');
paths.curation = fullfile(paths.workflowInputs, 'curation');
paths.curationReferenceModel = fullfile(paths.curation, 'shared', 'ThermoBuolVitCarbonConstrained.mat');
if ~isfile(paths.curationReferenceModel) && ~isempty(paths.archiveRoot)
    paths.curationReferenceModel = fullfile(paths.archiveWorkflowInputs, 'curation', 'shared', 'ThermoBuolVitCarbonConstrained.mat');
end
paths.dimes = fullfile(paths.workflowInputs, 'analyses', 'dimes');
paths.imm = fullfile(paths.workflowInputs, 'analyses', 'imm');
paths.nis = fullfile(paths.workflowInputs, 'analyses', 'nis');
paths.reconstruction = fullfile(paths.workflowInputs, 'reconstruction');
paths.reconstructionSequences = paths.reconstruction;
if ~isfolder(fullfile(paths.reconstructionSequences, 'FAAs')) && ~isempty(paths.archiveRoot)
    paths.reconstructionSequences = fullfile(paths.archiveWorkflowInputs, 'reconstruction');
end
paths.experimental = '';
if ~isempty(paths.archiveWorkflowInputs)
    paths.experimental = fullfile(paths.archiveWorkflowInputs, 'experimental');
end

paths.models.ptygCuratedTFA = fullfile(repositoryRoot, 'models', 'final', ...
    'curated_exchanges', 'PTYG');
paths.models.m9CuratedTFA = fullfile(repositoryRoot, 'models', 'final', ...
    'curated_exchanges', 'M9');
paths.models.ptygFullTFA = fullfile(repositoryRoot, 'models', 'final', ...
    'full_exchanges', 'PTYG');
paths.models.m9FullTFA = fullfile(repositoryRoot, 'models', 'final', ...
    'full_exchanges', 'M9');

% Current DiMEs and NIS scripts use the curated-exchange collections.
paths.models.richTFA = paths.models.ptygCuratedTFA;
paths.models.m9TFA = paths.models.m9CuratedTFA;

paths.external.matTFA = getenv('MATTFA_ROOT');
paths.external.ReMIND = getenv('REMIND_ROOT');
paths.external.redGEM = getenv('REDGEM_ROOT');
paths.external.NIS = getenv('NIS_ROOT');
paths.external.COBRA = getenv('COBRA_ROOT');
paths.external.RAVEN = getenv('RAVEN_ROOT');
paths.external.CPLEX = getenv('CPLEX_ROOT');
paths.external.Gurobi = getenv('GUROBI_ROOT');
end
