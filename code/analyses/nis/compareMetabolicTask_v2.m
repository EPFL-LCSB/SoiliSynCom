function compareMetabolicTask(output_PATH,sameBiomass,jaccard,names4biomass)
if nargin < 2
    sameBiomass = 1; % 1 if all models have the same biomass equation; 0 else
end

if nargin < 3
    jaccard = 1; % 1 for jaccard coefficient; 0 for %overlap
end

%% colormap for plotting
map_blue = [0.968627450980392,0.968627450980392,0.968627450980392;0.819607843137255,0.898039215686275,0.941176470588235;0.572549019607843,0.772549019607843,0.870588235294118;0.262745098039216,0.576470588235294,0.764705882352941;0.129411764705882,0.400000000000000,0.674509803921569;0.0196078431372549,0.188235294117647,0.380392156862745];
map_blue =imresize(map_blue,[100,3]);
map =imresize(redbluecmap,[30,3]);
map = min(max(map, 0), 1);

map_red = [0.968627450980392,0.968627450980392,0.968627450980392;0.992156862745098,0.858823529411765,0.780392156862745;0.956862745098039,0.647058823529412,0.509803921568627;0.839215686274510,0.376470588235294,0.301960784313725;0.698039215686275,0.0941176470588235,0.168627450980392;0.403921568627451,0,0.121568627450980];
map_red =imresize(map_red,[64,3]);
map_red = min(max(map_red, 0), 1);
map_combined = [flipud(map_blue); map_red]; % Flip blue to match red's gradient

%%
cd(output_PATH)

%% load data
folderNames = dir;
for i = 1:length(folderNames)
    names{i,1} = folderNames(i).name;
end
% get only relevant results
names = names(contains(names,'D0_Biosynthesis_PHB'));
% get results always on the same order
names = unique(names);
% id = names;
for i = 1:length(names)
    aux = regexp(names{i},'(.*)_D0_Biosynthesis_(.*)$','tokens');
    id{i,1} = aux{1}{1};
    Tasks{i,1} = aux{1}{2};
end
Tasks = unique(Tasks);
id = unique(id);

all_c = {};
all_alts = [];
all_sizeAlt = [];
meanSize = [];
%Tasks = Tasks(contains(Tasks,'polyP'));
for t = 1:length(Tasks)

    for i = 1:length(id)
        myFolder = strcat('./',id{i},'_D0_Biosynthesis_',Tasks{t},'/Models/',id{i},'/');
        filePattern = fullfile(myFolder, '*.mat');
        matFiles = dir(filePattern);
        if ~isempty(matFiles)
            if length(matFiles)>1

                var = strrep(matFiles(1).name,'.mat','');
                info = load(strcat(matFiles(1).folder,'/',matFiles(1).name));
            else
                var = strrep(matFiles.name,'.mat','');
                info = load(strcat(matFiles.folder,'/',matFiles.name));
            end
            info = info.(var);
            info = info.info_LMDP;
            %assignin('base',strcat('lump_info',num2str(i)),info);
            lump_info{i,t} = info;

            myFolder = strcat('./',id{i},'_D0_Biosynthesis_',Tasks{t},'/WorkSpaces/',id{i},'/',id{i},'/');
            filePattern = fullfile(myFolder, '*addBIOMASSnis_6.mat');
            if length(matFiles)>1

                var = strrep(matFiles(1).name,'.mat','');
                moreInfo = load(strcat(matFiles(1).folder,'/',matFiles(1).name),'LumpCstModel','DPsAll','bbb_metnames','summaryOfLUMPS');
            else
                matFiles = dir(filePattern);
                var = strrep(matFiles.name,'.mat','');
                moreInfo = load(strcat(matFiles.folder,'/',matFiles.name),'LumpCstModel','DPsAll','bbb_metnames','summaryOfLUMPS');
            end
            moreInfo_lump{i,t} = moreInfo;
        else
            lump_info{i,t} = 'NA';
            moreInfo_lump{i,t} ='NA';
        end
    end
end

%% get statistics

%U = {};
all_alts = [];
all_sizeAlt = [];
meanSize = [];

for t = 1:length(Tasks)
    sizeAlt = [];
    alts = [];

    for idx = 1:length(id)
        info1 = lump_info{idx,t};
        if iscell(info1)
            alts(idx,1) = size(info1,1);
            aux = [];
            for i=1:alts(idx,1)
                aux(i) = size(info1{i,4},1);
            end
            sizeAlt(idx,1) = min(aux);
            aux = info1(:,4);
            %U{idx,t} = length(unique(vertcat(aux{:})));
        else
            alts(idx,1) = 0;
            sizeAlt(idx,1) = 0;
            %U{idx,idx2} = 0;

        end
    end
    %all_c = [all_c  U];
    all_alts = [all_alts alts];
    all_sizeAlt = [all_sizeAlt sizeAlt];
end
meanSize = mean(all_sizeAlt,2);



% plot
yvalues = Tasks;
yvalues = strrep(yvalues,'_','-');
[~,idx] = ismember(order,id);
%yvalues = info.metNames(idx); %to get the names instead of ids
figure(1)
id = strrep(id,'_','-');
h = heatmap(id(idx),yvalues,all_alts(idx,:)','FontName','Avenir Next');
h.Title = 'Alternatives';
h.XLabel = 'Organisms';
%h.YLabel = 'BBBs';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_blue;
caxis([0,max(all_alts(:))])
saveas(h,'Alternatives','svg')
T = array2table(all_alts', 'VariableNames', id, 'RowNames', yvalues);
filename = strcat('statistics_Tasks_numAlt.xlsx');
writetable(T, filename, 'WriteRowNames', true);

figure(2)
h = heatmap(id(idx),yvalues,all_sizeAlt(idx,:)','FontName','Avenir Next');
h.Title = 'Subnetwork size';
h.XLabel = 'Organisms';
%h.YLabel = 'BBBs';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_red;
caxis([0,max(all_sizeAlt(:))])
saveas(h,'Subnetwork_size','svg')

T = array2table(all_sizeAlt', 'VariableNames', id, 'RowNames', yvalues);
filename = strcat('statistics_Tasks_sizeSub.xlsx');
writetable(T, filename, 'WriteRowNames', true);




figure(3)
h = heatmap('glycogen',id(idx),meanSize(idx),'FontName','Avenir Next');
h.Title = 'Subnetwork mean size';
h.YLabel = 'Organisms';
%h.YLabel = 'BBB category';
s = struct(h);
%s.XAxis.TickLabelRotation = 90;
h.Colormap = map_red;
caxis([0,max(meanSize(:))])
saveas(h,'Subnetwork_mean_size','svg')

T = array2table(meanSize, 'VariableNames', {'glycogen'}, 'RowNames', id);
filename = strcat('statistics_biosynthesis_MeanSize.xlsx');
writetable(T, filename, 'WriteRowNames', true);
%% comparison

for k = 1:length(Tasks)

    for idx = 1:length(id)
        for idx2 = 1:length(id)
            info1 = lump_info{idx,k};
            info2 = lump_info{idx2,k};


            if iscell(info1) && iscell(info2)
                aux1 = info1(:,4);
                aux1 = unique(vertcat(aux1{:}));

                aux2 = info2(:,4);
                aux2 = unique(vertcat(aux2{:}));

                overlap = numel(find(ismember(aux1,aux2)));
                union = numel(aux1) + numel(aux2) - overlap;
                if jaccard
                    J_lumps(idx,idx2) = overlap./union;
                else
                    J_lumps(idx,idx2) = overlap./numel(aux1) ;
                end
                clear overlap union aux1 aux2
            else
                J_lumps(idx,idx2) = 0;
            end

        end

        J_Tasks{k} = J_lumps;
    end

end



%save all data to a D0_Biosynthesis file
save('comparison_bioynthesis.mat','J_Tasks', 'id','Tasks')
id = strrep(id,'_','-');
Tasks = strrep(Tasks,'_','-');
%save average in xlsx for visualization and already plot a heatmap
for i = 1:length(Tasks)
    aux = J_Tasks{i};
    if ~isempty(aux)


        cl = clustergram(aux);
        r = cell2dbl(cl.RowLabels);

        figure(3+i)
        h = heatmap(id(r),id(r),aux(r,r),'FontName','Avenir Next');
        h.Title = Tasks{i};
        h.XLabel = 'Organisms';
        h.YLabel = 'Organisms';
        s = struct(h);
        s.XAxis.TickLabelRotation = 90;
        h.Colormap = map_red;
        saveas(h,strcat('similairity_',Tasks{i}),'svg')

        T = array2table(aux, 'VariableNames', id, 'RowNames', id);
        filename = strcat('similarity_',Tasks{i},'.xlsx');
        writetable(T, filename, 'WriteRowNames', true);
    end
end


%% extract the cost of biosynthesis for each one of the biomass building blocks for each organism

for idx = 1: length(id)
    cost = {};
    model_all = createModel();

    for i = 1 : length(Tasks)
        info1 = lump_info{idx,i};
        if iscell(info1)
        more = moreInfo_lump{idx,i};

        rxns1 = info1(:,1);
        formulas1 = info1(:,3);
        model1 = createModel();

        %in case we have metabolites in different compartments
        %for example in the yeast models

        if find(contains(formulas1,'_m'))
            m(i,1) = 1;
            formulas1 = strrep(formulas1,'_c','');
            formulas1 = strrep(formulas1,'_m','');
        end
        for k = 1:length(rxns1)
            model1 = addReactionTFA(model1,rxns1{k},formulas1{k});
        end
        aux = strsplit(Tasks{i},'-');
        aux = aux{1};
        aux = strcat(aux,'_c');
        coeff = [];
            f = find(ismember(model1.mets,aux));
            if f
            coeff = full(model1.S(f,:));
            model1.S = model1.S./coeff;
            end
        if not(ismember(info1(1,2),model1.mets)) && isfield(more,'DPsAll')
            more1 = moreInfo_lump{idx};
            f = find(ismember(more1.LumpCstModel.mets,info1(1,2)));
            idxNo = [find(more1.LumpCstModel.isTrans);  find(contains(more1.LumpCstModel.rxns,'EX_'))]; %need to fix this
            % get drains and transporters
            S = find(more1.LumpCstModel.S(f,:));
            S = setdiff(S,idxNo);
            f = find(ismember(more1.bbb_metnames,info1(1,2)));
            DPs = more1.DPsAll{f};
            indF = find(ismember(more1.LumpCstModel.varNames,strcat('F_',more1.LumpCstModel.rxns(S))));
            indR = find(ismember(more1.LumpCstModel.varNames,strcat('R_',more1.LumpCstModel.rxns(S))));
            DPs = DPs(indF,:) - DPs(indR,:);
            for d = 1:size(DPs,2)
                f = find(DPs(:,d));
            end
            f = unique(f);
            index = S(f);
            rxns = more1.LumpCstModel.rxns(index);
            %maybe it is not a transporter that is missing and the
            %final reaction is in the core already
            if ~isempty(rxns)
                formulas = printRxnFormula(more1.LumpCstModel,rxns);
                for k = 1:length(rxns)
                    model1 = addReactionTFA(model1,rxns{k},formulas{k});
                end
                coeff1 = sum(model1.S(:,end-length(rxns)+1:end),2);
                model1.S(:,end+1) = coeff1;
                model1.rxns(end+1) = strcat('aux_',info1(1,2)); model1.lb(end+1) = -100; model1.ub(end+1) = 100;
                if model1.S(find(ismember(model1.mets,info1(1,2))),end)<0
                    model1.S(:,end) = -1*model1.S(:,end);
                end
                model1 = removeRxns(model1,rxns);
            end
        end
        coeff1 = sum(model1.S,2);
        coeff1 = coeff1/length(rxns1);
        model1.S(:,end+1) = coeff1;
        model1.rxns(end+1) = info1(1,2);
        cost(i,1) = printRxnFormulaCost(model1,model1.rxns(end),0);

        else
            cost(i,1) = {'NA'};
        end
        %         else
        %             cost(i,1) = {':'};
        %         end
    end
    all_costs{idx,1} = cost;
%     T = array2table(cost, 'RowNames', Tasks);
%     filename = strcat('cost_',id{idx},'.xlsx');
%     writetable(T, filename, 'WriteRowNames', true);
    for t = 1:length(Tasks)
        model_all = addReactionTFAnis(model_all,strcat(id{idx},'_',Tasks{t}),cost{t,1});
    end
    all_rxns{idx} = model_all;
end
%save('./all_costs.mat','all_costs','all_rxns')
%% check if we have the same BBBs in the average lump
for t = 1:length(Tasks)
    J_cost = [];
    for idx = 1:length(id)
        model1 = all_rxns{idx};
        f1 = find(ismember(model1.rxns,strcat(id{idx},'_',Tasks(t))));
        m1 = find(model1.S(:,f1));

        for idx2 = 1:length(id)
            model2 = all_rxns{idx2};
            f2 = find(ismember(model2.rxns,strcat(id{idx2},'_',Tasks(t))));

            m2 = find(model2.S(:,f2));

            overlap = numel(find(ismember(m1,m2)));
            union = numel(m1) + numel(m2) - overlap;
            if jaccard
                J_cost(idx,idx2) = overlap./union;
            else
                J_cost(idx,idx2) = overlap./numel(m1) ;
            end
        end
        J_Tasks_cost{t,1} = J_cost;
    end
end

save('comparison_bioynthesis_cost.mat','J_Tasks_cost', 'id','Tasks')
id = strrep(id,'_','-');
Tasks = strrep(Tasks,'_','-');
%save average in xlsx for visualization and already plot a heatmap
for i = 1:length(Tasks)
    aux = J_Tasks_cost{i,1};
    if ~isempty(aux)


        cl = clustergram(aux);
        r = cell2dbl(cl.RowLabels);

        figure()
        h = heatmap(id(r),id(r),aux(r,r),'FontName','Avenir Next');
        h.Title = Tasks{i};
        h.XLabel = 'Organisms';
        h.YLabel = 'Organisms';
        s = struct(h);
        s.XAxis.TickLabelRotation = 90;
        h.Colormap = map_red;
        saveas(h,strcat('cost_similairity_',Tasks{i}),'svg')

        T = array2table(aux, 'VariableNames', id, 'RowNames', id);
        filename = strcat('cost_similarity_',Tasks{i},'.xlsx');
        writetable(T, filename, 'WriteRowNames', true);
    end
end
%%  now just get the coeff for ATP, NADH, NADPH, NAD and NADP
aux = [];
mets = {'PHB_c';'PHV_c';'cpd00002_c';'cpd00008_c';'C00003_c';'C00004_c';'C00006_c';'C00005_c'};
metNames = {'PHB';'PHV';'ATP';'ADP';'NAD';'NADH' ;'NADP';'NADPH'};
for t = 1:length(Tasks)
    coeff = [];
    for idx = 1:length(id)
        model1 = all_rxns{idx};
        f1 = find(ismember(model1.rxns,strcat(id{idx},'_',Tasks(t))));
        [~,fm] = ismember(mets,model1.mets);
        for i = 1:length(mets)
            if fm(i)
                aux(i,1) = full(model1.S(fm(i),f1));
            else
                aux(i,1) = 0;
            end
        end
        coeff= [coeff ; aux'];
  
    end
          Coeff_Tasks{t,1} = coeff;
end

save('comparison_bioynthesis_costCoeff.mat','Coeff_Tasks', 'mets','id','metNames')
id = strrep(id,'_','-');
Tasks = strrep(Tasks,'_','-');
%save average in xlsx for visualization and already plot a heatmap
for i = 1:length(Tasks)
    aux = Coeff_Tasks{i,1};
    if ~isempty(aux)


        cl = clustergram(aux);
        r = cell2dbl(cl.RowLabels);

        figure()
        aux = -aux;
        for k = 1:size(aux,2)

            f = find(aux(:,k) == -0);
            aux(f,k) = 0;
        end
        h = heatmap(metNames,id(r),aux(r,:),'FontName','Avenir Next');
        h.Title = Tasks{i};
        %h.XLabel = 'Organisms';
        h.YLabel = 'Organisms';
        s = struct(h);
        s.XAxis.TickLabelRotation = 90;
        h.Colormap = map_red;
        %saveas(h,strcat('costCoeff_',Tasks{i}),'svg')
% 
%         T = array2table(-aux, 'VariableNames', metNames, 'RowNames', id);
%         filename = strcat('costCoeff_',Tasks{i},'.xlsx');
%         writetable(T, filename, 'WriteRowNames', true);
    end
end
