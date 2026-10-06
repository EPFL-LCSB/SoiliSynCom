function reconstruct_draft_models(inputRoot, outputRoot, keggOutputRoot)
%RECONSTRUCT_DRAFT_MODELS Build the 21 SoilSynCom draft GEMs with RAVEN.
%
% inputRoot contains the bundled FAAs/. keggOutputRoot contains the retained
% KEGG search-output folders (BBET, BUOL, ...). The e-value threshold and
% scoring arguments reproduce the retained reconstruction commands. Each MAT
% file stores a variable named by its strain identifier.

paths = soilSynComPaths();

if nargin < 1 || isempty(inputRoot)
    inputRoot = paths.reconstructionSequences;
end
if nargin < 2 || isempty(outputRoot)
    outputRoot = fullfile(paths.outputRoot, 'draft_models');
end
if nargin < 3 || isempty(keggOutputRoot)
    keggOutputRoot = getenv('SOILSYNCOM_KEGG_OUTPUT_ROOT');
end
if isempty(keggOutputRoot) && ~isempty(paths.archiveRoot)
    keggOutputRoot = fullfile(paths.archiveRoot, 'reconstruction', 'KEGG_OUTPUTs');
end
if isempty(keggOutputRoot) || ~isfolder(keggOutputRoot)
    error(['Provide keggOutputRoot or set SOILSYNCOM_KEGG_OUTPUT_ROOT to ' ...
        'the extracted KEGG_OUTPUTs directory.']);
end
if ~isfolder(outputRoot)
    mkdir(outputRoot);
end

models = {
    'bbet', 'br_betae.fasta',                'BBET';
    'buol', 'burkholderia_olga.fasta',       'BUOL';
    'cauj', 'caulobacter_ji38.fasta',        'CAUJ';
    'cpi',  'ch_pinensis.fasta',             'CPI';
    'cohm', 'cohnella_mfer.fasta',           'COHM';
    'cuvk', 'curtobacterium_vkm.fasta',      'CUVK';
    'cxy',  'c_xylanilytica.fasta',          'CXY';
    'dri',  'de_riboflavina.fasta',          'DRI';
    'flk',  'flavobacterium_kbs0721.fasta',  'FLK';
    'lur',  'lu_rhozovicinus.fasta',         'LUR';
    'lgu',  'ly_gummosus.fasta',             'LGU';
    'mlot', 'mesorhizobium_loti.fasta',      'MLOT';
    'mgk',  'm_ginsenosidivorax_khi28.fasta','MGK';
    'mip',  'microbacterium_pamc.fasta',     'MIP';
    'pzu',  'ph_zucineum_hlk1.fasta',        'PZU';
    'pscf', 'pseudomonas_cfsan.fasta',       'PSCF';
    'pmor', 'ps_moraviensis.fasta',          'PMOR';
    'racd', 'rahnella_cdc.fasta',            'RACD';
    'rfa',  'rh_fascians.fasta',             'RFA';
    'tro',  'ta_robiniae.fasta',             'TRO';
    'vpe',  'v_paradoxus_eps.fasta',         'VPE'
};

for index = 1:size(models, 1)
    strain = models{index, 1};
    fasta = fullfile(inputRoot, 'FAAs', models{index, 2});
    keggOutput = fullfile(keggOutputRoot, models{index, 3}, 'kegg_output');
    draft = getKEGGModelForOrganism( ...
        strain, fasta, 'dataDir', keggOutput, 0, 0, 0, 1e-20, 0.8, 0.3, -1);
    payload = struct();
    payload.(strain) = draft;
    save(fullfile(outputRoot, [strain '.mat']), '-struct', 'payload');
end
end
