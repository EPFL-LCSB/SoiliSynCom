function getConnections(output_PATH,D)
cd(output_PATH)

%% load data
folderNames = dir;
for i = 1:length(folderNames)
    names{i,1} = folderNames(i).name;
end
% get results for a specific D
names = names(contains(names,strcat('_D',num2str(D),'_core')));
% get results always on the same order
names = unique(names);
id = names;
id = strrep(names,strcat('_D',num2str(D),'_core'),'');


for i = 1:length(names)
    myFolder = strcat('./',names{i},'/WorkSpaces/',id{i},'/',id{i},'/');
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    info = load(strcat(matFiles.folder,'/',matFiles.name),'selectedPaths','GSM_ForAdjMat','core_ss');
    eval(strcat('info',num2str(i),'= info'));
end

%% get connections between subsystems
for idx = 1:length(names)
    info = eval(strcat('info',num2str(idx)));
a = [];
a1 = [];
for i = 1:length(info.core_ss)
    
    for j = 1:length(info.core_ss)
        
        a(:,1) = info.selectedPaths.joiningReacs{i,j,D};
        a1(:,1) = info.selectedPaths.joiningReacs{j,i,D};
        
        a = [a;a1];
        a = unique(a);
        a = info.GSM_ForAdjMat.rxns(a);
        num_connecting_reactions(i,j) = length(a);
        
        b = concatenateList(a,'|');
        NIS_connecting_reactions{i,j} = b;
        
        clear a a1
        
    end
    
    clear a a1
    
end
T = array2table(num_connecting_reactions, 'VariableNames', info.core_ss, 'RowNames', info.core_ss);
filename = strcat('NIS_connecting_reactions_statistics_',id{idx},'.xlsx');
writetable(T, filename, 'WriteRowNames', true);


T = array2table(NIS_connecting_reactions, 'VariableNames', info.core_ss, 'RowNames', info.core_ss);
filename = strcat('NIS_connecting_reactions_',id{idx},'.xlsx');
writetable(T, filename, 'WriteRowNames', true);


a = [];
a1 = [];
for i = 1:length(info.core_ss)
    
    for j = 1:length(info.core_ss)
        
        a(:,1) = info.selectedPaths.joiningMets{i,j,D};
        a1(:,1) = info.selectedPaths.joiningMets{j,i,D};
        
        a = [a;a1];
        a = unique(a);
        a = info.GSM_ForAdjMat.mets(a);
        num_connecting_mets(i,j) = length(a);
        
        b = concatenateList(a,'|');
        NIS_connecting_metabolites{i,j} = b;
        
        clear a a1
        
    end
    
    clear a a1
    
end

T = array2table(NIS_connecting_metabolites, 'VariableNames', info.core_ss, 'RowNames', info.core_ss);
filename = strcat('NIS_connecting_metabolites_',id{idx},'.xlsx');
writetable(T, filename, 'WriteRowNames', true);

T = array2table(num_connecting_mets, 'VariableNames', info.core_ss, 'RowNames', info.core_ss);
filename = strcat('NIS_connecting_metabolites_statistics_',id{idx},'_D',num2str(D),'.xlsx');
writetable(T, filename, 'WriteRowNames', true);
end
end