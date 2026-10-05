function pathwayCompleteness_mets(output_PATH,data_path)
cd(output_PATH)
load(data_path)
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
D = 0;

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
    eval(strcat('info',num2str(i),'= info'));
end
clear info
%% calaculate the pathway completeness with respect to mets
 metabolites = metabolites([1,2,3,5]); %reorder subsystems
for j = 1:length(names)
    for i = 1:length(info1.initialSubSystems)-1
        aux = eval(strcat('info',num2str(j)));
        mets = aux.mets_ss_ID{i};
        mets = strrep(mets,'_c','');
        mets = strrep(mets,'_e','');
        mets = strrep(mets,'_p','');
        mets = unique(mets);
        mets_ref = unique(metabolites{i});
        complete(i,j) = length(find(ismember(mets,mets_ref)))/length(mets_ref);
    end
end

% plot 
 
cl = clustergram(complete);
r_r = cell2dbl(cl.RowLabels);
r_c = cell2dbl(cl.ColumnLabels);
figure()
h_mets = heatmap(id(r_c),info1.initialSubSystems(r_r),complete(r_r,r_c),'FontName','Avenir Next');
h_mets.Title = 'Pathway completeness'; h_mets.XLabel = 'Organisms';
h_mets.YLabel = 'Pathways'; s = struct(h_mets); 
s.XAxis.TickLabelRotation= 90; h_mets.Colormap = map_red; h_mets.CellLabelFormat = '%.2f';
caxis([0,1])
saveas(h_mets,'pth_complete','svg') 
T =array2table(complete, 'VariableNames', id, 'RowNames',info1.initialSubSystems(1:4)); 
filename =strcat(output_PATH,'/pth_complete.xlsx');
writetable(T, filename, 'WriteRowNames', true);
end
