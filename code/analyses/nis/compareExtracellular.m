function compareExtracellular(output_PATH,onlySecretions,media)
if nargin<2
   onlySecretions = 1;
end
% you can exclude the media components from the comparison
if nargin<3
    media = {};
end
%this function asseses the similarity of the overall metabolite exchange
%subnetworks or the secretion profiles
cd(output_PATH)
%% colormaps
map_blue = [0.968627450980392,0.968627450980392,0.968627450980392;0.819607843137255,0.898039215686275,0.941176470588235;0.572549019607843,0.772549019607843,0.870588235294118;0.262745098039216,0.576470588235294,0.764705882352941;0.129411764705882,0.400000000000000,0.674509803921569;0.0196078431372549,0.188235294117647,0.380392156862745];
map_blue =imresize(map_blue,[100,3]);
map =imresize(redbluecmap,[30,3]);
map = min(max(map, 0), 1);
map_red = [0.968627450980392,0.968627450980392,0.968627450980392;0.992156862745098,0.858823529411765,0.780392156862745;0.956862745098039,0.647058823529412,0.509803921568627;0.839215686274510,0.376470588235294,0.301960784313725;0.698039215686275,0.0941176470588235,0.168627450980392;0.403921568627451,0,0.121568627450980];
map_red =imresize(map_red,[64,3]);
map_red = min(max(map_red, 0), 1);
%% load data
folderNames = dir;
for i = 1:length(folderNames)
    names{i,1} = folderNames(i).name;
end
names = names(contains(names,'_D0_EX'));
% get results always on the same order
names = unique(names);

id = names;
id = strrep(names,'_D0_EX','');


for i = 1:length(names)
    myFolder = strcat('./',names{i},'/Models/',id{i},'/');
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    info = load(strcat(matFiles.folder,'/',matFiles.name));
    var = strrep(matFiles.name,'.mat','');
    info = info.(var);
    eval(strcat('info',num2str(i),'= info'));
end

%% get statistics
for idx = 1:length(names)
    for idx2 = 1:length(info1.info_Econnect)
        extra1 = eval(strcat('info',num2str(idx)));
        extra1= extra1.info_Econnect;
        if not(isempty(extra1{idx2,2}))
            alts(idx,idx2) = length(extra1{idx2,2});
            sizeAlt(idx,idx2) = length(extra1{idx2,2}{1});
            alts2(idx,idx2) = 1;
        else
            alts(idx,idx2) = 0;
            sizeAlt(idx,idx2) = 0;
            alts2(idx,idx2) = 0;
        end
    end
end

xvalues = id;
yvalues = info1.info_Econnect(:,1);


figure(1)
h = heatmap(id,yvalues,alts','FontName','Avenir Next');
h.Title = 'Alternatives';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_blue;
caxis([0,max(alts(:))])
T = array2table(alts', 'VariableNames', id, 'RowNames', yvalues);
filename = 'statistics_extracellular_a.xlsx';
writetable(T, filename, 'WriteRowNames', true);
figure(2)
h = heatmap(id,yvalues,sizeAlt','FontName','Avenir Next');
h.Title = 'Subnetwork size';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_red;
caxis([0,max(sizeAlt(:))])
T = array2table(sizeAlt', 'VariableNames', id, 'RowNames', yvalues);
filename = 'statistics_extracellular_b.xlsx';
writetable(T, filename, 'WriteRowNames', true);
%% get the subnetworks  for each model and carbon source
% first list all products in general
List = {};
for idx = 1:length(names)
    extra1 = eval(strcat('info',num2str(idx)));
    extra1 = extra1.info_Econnect;
    for idx2 = 1:length(info1.info_Econnect)
        aux1 = extra1{idx2,2};
        if not(isempty(aux1))
            aux1 = unique(vertcat(aux1{:}));
            if onlySecretions
            aux1 = aux1(contains(aux1,'EX_'));
            end
            aux1 = setdiff(aux1,[media;extra1{idx2,1}]);
            List = [List;aux1];
            List = unique(List);
        end
    end
end

%now check if the product (or reaction) is at least with one carbon source
for l=1:length(List)
for idx = 1:length(names)
    extra1 = eval(strcat('info',num2str(idx)));
    extra1 = extra1.info_Econnect;
    for idx2 = 1:length(info1.info_Econnect)
        aux1 = extra1{idx2,2};
        if not(isempty(aux1))
            aux1 = unique(vertcat(aux1{:}));
            if onlySecretions
            aux1 = aux1(contains(aux1,'EX_'));
            end
            aux1 = setdiff(aux1,[media;extra1{idx2,1}]);
            if (ismember(List{l},aux1))
                isthere(l,idx) = 1;

            end
            end
        end
    end
end

% get the secretion data or reaction data
data = [];
i = 1;
for idx2 = 1:length(info1.info_Econnect)
    for idx = 1:length(names)
        extra1 = eval(strcat('info',num2str(idx)));
        extra1 = extra1.info_Econnect;
        aux1 = extra1{idx2,2};
        if not(isempty(aux1))
            aux1 = unique(vertcat(aux1{:}));
            if onlySecretions
            aux1 = aux1(contains(aux1,'EX_'));
            end
            aux1 = setdiff(aux1,[media;extra1{idx2,1}]);
            data = [data ismember(List,aux1)];
        else
            data = [data zeros(length(List),1)];
        end
        i=i+1;
    end
    
end
id_aux = extra1(:,1);
id_aux = strrep(id_aux,'EX_','');
id_aux = strrep(id_aux,'_e','');
new_tags = {};
for i = 1:length(id_aux)
    new_tags = [new_tags ;strcat(id_aux{i},'-',id)];
end
List2 = strrep(List,'EX_','');
% save secretion profiles
T = array2table(data, 'VariableNames', new_tags, 'RowNames', List2);
filename = 'subnetwork_profiles.xlsx';
writetable(T, filename, 'WriteRowNames', true);

% similarity of secretion profiles and PCo 
s = sum(data,1);
aux = [1:length(data)];
aux = setdiff(aux,find(s==0));
data3 = data(:,aux);
similarityMatrix = computeJaccardSimilarityMatrix(data3);
% save secretion data ready for PCo
T = array2table(similarityMatrix, 'VariableNames', new_tags(aux), 'RowNames', new_tags(aux));
filename = 'subnetwork_similarity.xlsx';
writetable(T, filename, 'WriteRowNames', true);
