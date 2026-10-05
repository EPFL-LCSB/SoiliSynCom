function compareCore(output_PATH,D,jaccard)
if nargin < 3
    jaccard = 1; % 1 for jaccard coefficient; 0 for %overlap
end
cd(output_PATH)
%% colormap for plotting
map_red = [0.968627450980392,0.968627450980392,0.968627450980392;0.992156862745098,0.858823529411765,0.780392156862745;0.956862745098039,0.647058823529412,0.509803921568627;0.839215686274510,0.376470588235294,0.301960784313725;0.698039215686275,0.0941176470588235,0.168627450980392;0.403921568627451,0,0.121568627450980];
map_red =imresize(map_red,[64,3]);
map_red = min(max(map_red, 0), 1);
map_blue = [0.968627450980392,0.968627450980392,0.968627450980392;0.819607843137255,0.898039215686275,0.941176470588235;0.572549019607843,0.772549019607843,0.870588235294118;0.262745098039216,0.576470588235294,0.764705882352941;0.129411764705882,0.400000000000000,0.674509803921569;0.0196078431372549,0.188235294117647,0.380392156862745];
map_blue =imresize(map_blue,[64,3]);
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
    myFolder = strcat('./',names{i},'/Models/',id{i},'/');
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    info = load(strcat(matFiles.folder,'/',matFiles.name));
    var = strrep(matFiles.name,'.mat','');
    info = info.(var);
    WS=fullfile(strcat('./',names{i},'/WorkSpaces/',strrep(names{i},'_D0_core',''),...
    '/',strrep(names{i},'_D0_core',''),'/'),'*.mat');
    WSMat=dir(WS);
    SS = load(strcat(WSMat.folder,'/',WSMat.name),'rxns_ss','mets_ss','GSM_ForAdjMat');
    rxns_coreall=SS.rxns_ss;
    mets_coreall=SS.mets_ss;
    for v=1:length(rxns_coreall)
        rxns_SSall{v,1} = SS.GSM_ForAdjMat.rxns(rxns_coreall{v});
    end
    for v=1:length(mets_coreall)
        mets_SSall{v,1} = SS.GSM_ForAdjMat.mets(mets_coreall{v});
    end

    info.mets_ss_ID=mets_SSall;
    info.rxns_ss_ID=rxns_SSall;
    eval(strcat('info',num2str(i),'= info'));
end
clear info
%% initial pathways
% number of reactions and metabolites
if D == 0
for j = 1:length(names)
    for i = 1:length(info1.initialSubSystems)
        aux = eval(strcat('info',num2str(j)));
        mets(j,i) = length(aux.mets_ss_ID{i});
        rxns(j,i) = length(aux.rxns_ss_ID{i});
    end
end

% plot
figure(1)
h_mets = heatmap(id,info1.initialSubSystems,mets','FontName','Avenir Next');
h_mets.Title = 'Number of Metabolites';
h_mets.XLabel = 'Organisms';
h_mets.YLabel = 'Pathways';
s = struct(h_mets);
s.XAxis.TickLabelRotation = 90;
h_mets.Colormap = map_red;
h_mets.CellLabelFormat = '%.2f';
caxis([0,max(mets(:))])
T = array2table(mets', 'VariableNames', id, 'RowNames', info1.initialSubSystems);
filename = 'statistics_mets_initial_subsystems.xlsx';
writetable(T, filename, 'WriteRowNames', true);

figure(2)
h_rxns = heatmap(id,info1.initialSubSystems,rxns','FontName','Avenir Next');
h_rxns.Title = 'Number of Reactions';
h_rxns.XLabel = 'Organisms';
h_rxns.YLabel = 'Pathways';
s = struct(h_rxns);
s.XAxis.TickLabelRotation = 90;
h_rxns.Colormap = map_blue;
h_rxns.CellLabelFormat = '%.2f';
caxis([0,max(rxns(:))])
T = array2table(rxns', 'VariableNames', id, 'RowNames', info1.initialSubSystems);
filename = 'statistics_rxns_initial_subsystems.xlsx';
writetable(T, filename, 'WriteRowNames', true);

%% similarity pathways
ptrn = '_[a-z]+';
L = length(info1.initialSubSystems);
for l = 1:L
for idx = 1:length(names)
    for idx2 = 1:length(names)
        % metabolites
        path1 = eval(strcat('info',num2str(idx)));
        path2 = eval(strcat('info',num2str(idx2)));        
            a = path1.mets_ss_ID{l};
            l_aa = length(a);
            b = path2.mets_ss_ID{l};
            overlap_M = numel(find(ismember(a,b)));
            union_M = numel(a) + numel(b) - overlap_M;

            a = regexprep(a,ptrn,'');
            a = unique(a);
            l_aa_n = length(a);
            b = regexprep(b,ptrn,'');
            b = unique(b);

            overlap_M_n = numel(find(ismember(a,b)));
            union_M_n = numel(a) + numel(b) - overlap_M_n;

        %reactions
            a = path1.rxns_ss_ID{l};
            l_a = length(a);
            b = path2.rxns_ss_ID{l};
            overlap = numel(find(ismember(a,b)));
            union = numel(a) + numel(b) - overlap;
       

        if  jaccard
            J_core_rxns(idx,idx2) = overlap./union;
        else
            J_core_rxns(idx,idx2) = overlap./l_a;
        end
if l==L
        % here we find all the metabolites involed in ETC
        rxns1 = a;
        [~,idx_r] = ismember(rxns1,path1.OriginalGEM.rxns);
        a = {};
        for i = 1:length(idx_r)
            a = [a; path1.OriginalGEM.mets(find(path1.OriginalGEM.S(:,idx_r(i))))];
        end
        a = unique(a);

        rxns2 = b;
        [~,idx_r] = ismember(rxns2,path2.OriginalGEM.rxns);
        b = {};
        for i = 1:length(idx_r)
            b = [b; path2.OriginalGEM.mets(find(path2.OriginalGEM.S(:,idx_r(i))))];
        end
        b = unique(b);

        overlap_M = numel(find(ismember(a,b)));
        union_M = numel(a) + numel(b) - overlap_M;
        l_aa =  length(unique(a));

        a = regexprep(a,ptrn,'');
        a = unique(a);
        l_aa_n =  length(unique(a));

        b = regexprep(b,ptrn,'');
        b = unique(b);

        overlap_M_n = numel(find(ismember(a,b)));
        union_M_n = numel(a) + numel(b) - overlap_M_n;
end
        if  jaccard
            J_core_mets(idx,idx2) = overlap_M./union_M;
            J_core_mets_n(idx,idx2) = overlap_M_n./union_M_n;
        else
            J_core_mets(idx,idx2) = overlap_M./l_aa;
            J_core_mets_n(idx,idx2) = overlap_M_n./l_aa_n;
        end


    end
end
info1.initialSubSystems = strrep(info1.initialSubSystems,'/','_');
T = array2table(J_core_rxns, 'VariableNames', id, 'RowNames', id);
filename = strcat(info1.initialSubSystems{l},'_rxns_similarity.xlsx');
writetable(T, filename, 'WriteRowNames', true);
T = array2table(J_core_mets_n, 'VariableNames', id, 'RowNames', id);
filename = strcat(info1.initialSubSystems{l},'_mets_similarity.xlsx');
writetable(T, filename, 'WriteRowNames', true);
figure()
h_mets = heatmap(id,id,J_core_mets_n,'FontName','Avenir Next');
h_mets.Title = strcat('Metabolite Similarity - ',info1.initialSubSystems{l});
h_mets.CellLabelFormat = '%.2f';
s = struct(h_mets);
s.XAxis.TickLabelRotation = 90;
h_mets.Colormap = map_red;
caxis([0,1])

figure()
h_rxns = heatmap(id,id,J_core_rxns,'FontName','Avenir Next');
h_rxns.CellLabelFormat = '%.2f';
h_rxns.Title = strcat('Reaction Similarity - ',info1.initialSubSystems{l});
s = struct(h_rxns);
s.XAxis.TickLabelRotation = 90;
h_rxns.Colormap = map_blue;
caxis([0,1])

end
end
%% for networks D0, D1, D2, etc
ptrn = '_[a-z]+';
for idx = 1:length(names)
    for idx2 = 1:length(names)
            clear overlap union
            D01 = eval(strcat('info',num2str(idx)));
            D02 = eval(strcat('info',num2str(idx2)));
            mets_ss_all1 = D01.core_model.mets;
            mets_ss_all2 = D02.core_model.mets;
            overlap = numel(find(ismember(mets_ss_all1,mets_ss_all2)));
            union = numel(mets_ss_all1) + numel(mets_ss_all2) - overlap;

            if jaccard 
            J_core_mets_all(idx,idx2) = overlap./union;
            else
            J_core_mets_all(idx,idx2) = overlap./numel(mets_ss_all1);
            end


            mets_ss_all1 = regexprep(mets_ss_all1,ptrn,'');
            mets_ss_all1= unique(mets_ss_all1);
            mets_ss_all2 = regexprep(mets_ss_all2,ptrn,'');
            mets_ss_all2= unique(mets_ss_all2);
            overlap = numel(find(ismember(mets_ss_all1,mets_ss_all2)));
            union = numel(mets_ss_all1) + numel(mets_ss_all2) - overlap;

            if jaccard 
            J_core_mets_all_n(idx,idx2) = overlap./union;
            else
            J_core_mets_all_n(idx,idx2) = overlap./numel(mets_ss_all1);
            end
            
            rxns_ss_all1 = D01.core_model.rxns;
            rxns_ss_all2 = D02.core_model.rxns;
            overlap = numel(find(ismember(rxns_ss_all1,rxns_ss_all2)));
            union = numel(rxns_ss_all1) + numel(rxns_ss_all2) - overlap;

            if jaccard 
            J_core_rxns_all(idx,idx2) = overlap./union;
            else
            J_core_rxns_all(idx,idx2) = overlap./numel(rxns_ss_all1);
            end
            
    end
end
T = array2table(J_core_rxns_all, 'VariableNames', id, 'RowNames', id);
filename = strcat('D',num2str(D),'_rxns_similarity.xlsx');
writetable(T, filename, 'WriteRowNames', true);
T = array2table(J_core_mets_all_n, 'VariableNames', id, 'RowNames', id);
filename = strcat('D',num2str(D),'_mets_similarity.xlsx');
writetable(T, filename, 'WriteRowNames', true);
figure()
h_mets = heatmap(id,id,J_core_mets_all_n,'FontName','Avenir Next');
h_mets.Title = strcat('Metabolite Similarity - ','D',num2str(D));
h_mets.CellLabelFormat = '%.2f';
s = struct(h_mets);
s.XAxis.TickLabelRotation = 90;
h_mets.Colormap = map_red;
caxis([0,1])

figure()
h_rxns = heatmap(id,id,J_core_rxns_all,'FontName','Avenir Next');
h_rxns.Title = strcat('Reaction Similarity - ','D',num2str(D));
h_rxns.CellLabelFormat = '%.2f';
s = struct(h_rxns);
s.XAxis.TickLabelRotation = 90;
h_rxns.Colormap = map_blue;
caxis([0,1])