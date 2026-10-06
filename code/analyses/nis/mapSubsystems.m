%% script to map reaction metabolic subsystems
paths = soilSynComPaths();
pathToGEMs = paths.models.richTFA;
outPath = fullfile(paths.outputRoot, 'NIS', 'SubSystems'); %in case you don't want to overwrite the models
% please provide different output path than the folder you have initialy the GEMs
% cd(outPath)
if ~isfolder(outPath)
    mkdir(outPath)
end
addSoilSynComDependencies(paths);
ids = 'KEGG'; %BiGG or KEGG for now
if strcmp(ids,'KEGG')
    load(fullfile(paths.nis, 'metabolites4subsystems_KEGG.mat'))
    
    met_pairs_to_remove = load(fullfile(paths.nis, 'met_pairs_to_remove_KEGG.mat'));
    met_pairs_to_remove = met_pairs_to_remove.met_pairs_to_remove;
    
    ss = {'Glycolysis/Gluconeogenesis';...
        'Pentose Phosphate Pathway';...
        'Citric Acid Cycle';
        'Pyruvate Metabolism'};
    
elseif strcmp(ids,'BiGG')
    error('Provide the BiGG subsystem-mapping MAT files before selecting ids = BiGG.');
end

% load models
% here you can read directly the models from the file or
cd(pathToGEMs)
folderNames = dir;
for i = 1:length(folderNames)
    names{i,1} = folderNames(i).name;
end
names = names(contains(names,'.mat'));
% get results always on the same order
names = unique(names);
% provide the model names
%names = {'model1';'model2'};

for idx_m = 1:length(names)
    data = load(names{idx_m}); % Load the .mat file into a struct
    modelName = fieldnames(data); % Get the variable name stored in the .mat file
    model = data.(modelName{1}); 
%     model=load(names{idx_m}).names{idx_m};
    %here we assume that the model structur is named 'model'
    if size(model.metCompSymbol,1)<size(model.metCompSymbol,2)
        model.metCompSymbol = model.metCompSymbol';
        if exist('ttmodel','var')
            ttmodel.metCompSymbol = model.metCompSymbol';
        end
    end
    model.subSystems = {};
    InorgMetSEEDIDs = findInorganicMets(model);
    
    %% map based on pairs
    
    model.subSystems = cell(length(model.rxns),1);
    model = checkTransportNIS(model);
    
    for i = 1:length(model.subSystems)
        model.subSystems{i,1} = '|';
    end
    
    %exclude exchange and transporters
    selExc = full((sum(abs(model.S)==1,1) ==1) & (sum(model.S~=0) == 1))';
    selExc = find(selExc);
    transport = find(model.isTrans);
    selExc = unique([selExc;transport]);
    
    for i = 1:size(ss,1)
        if i~=4 && i~=6 % if we don't want to map carbon fixation and glyoxylate metabolism
            selExci = selExc;
            modeltemp = model;
            met_pairs_to_removei = met_pairs_to_remove;
            
            if strcmp(ss{i},'Glycolysis/Gluconeogenesis')
                aux = {'C00074_c:C00022_c';'C00015_c:C00075_c';'C00044_c:C00035_c';...
                    'C00112_c:C00063_c';'pep_c:pyr_c';'udp_c:utp_c';'gtp_c:gdp_c';'cdp_c:ctp_c';...
                    'cpd00061_c:cpd00020_c';'cpd00014_c:cpd00062_c';'cpd00038_c:cpd00031_c';...
                    'cpd00096_c:cpd00052_c'};
                met_pairs_to_removei = setdiff(met_pairs_to_remove,aux);
                met_pairs_to_removei = [met_pairs_to_removei;{'dha_c:dhap_c';'C00184_c:C00111_c';'cpd00157_c:cpd00095_c'}];
                selExci = unique([selExci;find(ismember(model.rxns,findTransRxns(model)))]);
            end
            
            if strcmp(ss{i},'Pyruvate Metabolism')
                aux = {'C00074_c:C00022_c';'C00010_c:C00024_c';'C00015_c:C00075_c';...
                    'C00044_c:C00035_c';'C00112_c:C00063_c';'pep_c:pyr_c';'coa_c:accoa_c';...
                    'udp_c:utp_c';'gtp_c:gdp_c';'cdp_c:ctp_c';'cpd00061_c:cpd00020_c';...
                    'cpd00010_c:cpd00022_c';'cpd00014_c:cpd00062_c';'cpd00038_c:cpd00031_c';...
                    'cpd00096_c:cpd00052_c'};
                met_pairs_to_removei = setdiff(met_pairs_to_remove,aux);
                selExci = unique([selExci;find(ismember(model.rxns,findTransRxns(model)))]);
            end
            
            
            if strcmp(ss{i},'Citric Acid Cycle')
                aux = {'C00010_c:C00024_c';'C00010_c:C00091_c';'coa_c:accoa_c';...
                    'coa_c:succoa_c';'coa_m:accoa_m';'coa_m:succoa_m';...
                    'cpd00010_c:cpd00022_c';'cpd00010_c:cpd00078_c'};
                met_pairs_to_removei = setdiff(met_pairs_to_remove,aux);
                %met_pairs_to_removei = [met_pairs_to_removei;'ampnt_c:acampnt_c'];
                selExci = unique([selExci;find(ismember(model.rxns,findTransRxns(model)))]);
            end
            
            Smatrix = removeMetPairsFromS(model, met_pairs_to_removei, InorgMetSEEDIDs);
            modeltemp.S = Smatrix;
            modeltemp.mets = strrep(modeltemp.mets,strcat('_',modeltemp.metCompSymbol),'');
            
            for j = 1:length(model.rxns)
                if not(ismember(j,selExci))
                    mets = modeltemp.mets(find(modeltemp.S(:,j)));
                    if all(ismember(mets,metabolites{i})) && length(mets)>1
                        model.subSystems{j} = strcat(model.subSystems{j},ss{i},'|');
                    end
                end
            end
        end
    end
    f = find(ismember(model.subSystems,...
        {'|Glycolysis/Gluconeogenesis|Pentose Phosphate Pathway|';...
        '|Glycolysis/Gluconeogenesis|Pyruvate Metabolism|'}));
    model.subSystems(f) = {'Glycolysis/Gluconeogenesis'};
    f = find(ismember(model.subSystems,'|Citric Acid Cycle|Pyruvate Metabolism|'));
    model.subSystems(f) = {'Citric Acid Cycle'};
    f = find(ismember(model.rxns,'R00431'));
    model.subSystems(f) = {'Glycolysis/Gluconeogenesis'};
%     f = find(ismember(model.rxns,'R00024'));
%     model.subSystems(f) = {'Carbon Fixation'};
    
%     %for accumulibacter !!!
%     PHA = {'PHBphaB'
%         'PHBphaC'
%         'PHVphaA'
%         'PHVphaB'
%         'PHVphaC'
%         'PH2MVphaA'
%         'PH2MVphaB'
%         'PH2MVphaC'
%         'PH2MBphaA'
%         'PH2MBphaB'
%         'PH2MBphaC'
%         'PHBphaZ'
%         'PHVphaZ'
%         'PH2MVphaZ'
%         'PH2MBphaZ'};
%     f = find(ismember(model.rxns,PHA));
%     model.subSystems(f) = {'PHA'};
%     
%     polyP = {'PPK1'
%         'PPK2g'
%         'PPK2a'
%         'PPX'
%         'PAPg'
%         'PAPa'
%         'NADK'
%         'PPGK'};
%     
%     Glycogen = {'GS1_GlycI'
%         'GS2_GlycII'
%         'GS3_GlycIII'
%         'GD1_GlycIV'
%         'GD2_GlycV'
%         'R02111_mod'
%         'R02112_mod'
%         'R01790_mod'
%         'R00948'
%         'R11530'
%         'R01557'};
%     f = find(ismember(model.rxns,Glycogen));
%     model.subSystems(f) = {'Glycogen'};
%     
%     polyP = {'PPK1'
%         'PPK2g'
%         'PPK2a'
%         'PPX'
%         'PAPg'
%         'PAPa'
%         'NADK'
%         'PPGK'};
%     
%     f = find(ismember(model.rxns,polyP));
%     model.subSystems(f) = {'polyP'};
    
    model = generate_ETC_SubSyst(model, 'NA');
    model.subSystems = strrep(model.subSystems,'|','');
    if exist('ttmodel','var')
        ttmodel.subSystems = model.subSystems;
        save(strcat(outPath,'/',names{idx_m}),'ttmodel','model')
    else
        % Extract base filename without extension
        [~, varName, ~] = fileparts(names{idx_m});
        
        % Remove 'Thermo' and 'PTYGVitconstrained' from the name
        varName = erase(varName, {'Thermo', 'PTYGVitconstrained'});
        
        % Rename the model variable dynamically
        eval(['ttmodel_',varName ' = model;']);
        save(fullfile(outPath, names{idx_m}), strcat('ttmodel_',varName));
%         % Save using the cleaned variable name
%         save(strcat(outPath,'/',names{idx_m}),'model','model')
    end
end
