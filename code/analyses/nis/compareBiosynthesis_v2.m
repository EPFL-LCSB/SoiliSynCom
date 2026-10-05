function compareBiosynthesis(output_PATH,sameBiomass,jaccard,names4biomass)
if nargin < 2
    sameBiomass = 1; % 1 if all models have the same biomass equation; 0 else
end

if nargin < 3
    jaccard = 1; % 1 for jaccard coefficient; 0 for %overlap
end

%% colormap for plotting
map_blue = [0.968627450980392,0.968627450980392,0.968627450980392;0.819607843137255,0.898039215686275,0.941176470588235;0.572549019607843,0.772549019607843,0.870588235294118;0.262745098039216,0.576470588235294,0.764705882352941;0.129411764705882,0.400000000000000,0.674509803921569;0.0196078431372549,0.188235294117647,0.380392156862745];
map_blue =imresize(map_blue,[100,3]);
map_blue =imresize(map_blue,[30,3]);
map_blue = min(max(map_blue, 0), 1);

map_red = [0.968627450980392,0.968627450980392,0.968627450980392;0.992156862745098,0.858823529411765,0.780392156862745;0.956862745098039,0.647058823529412,0.509803921568627;0.839215686274510,0.376470588235294,0.301960784313725;0.698039215686275,0.0941176470588235,0.168627450980392;0.403921568627451,0,0.121568627450980];
map_red =imresize(map_red,[64,3]);
map_red = min(max(map_red, 0), 1);
%%
cd(output_PATH)

if sameBiomass
    names4biomass = {};
    idx_biomass = 0;
elseif  sameBiomass == 0 && nargin<4
    error('Models with different biomass compositions.Please define in the inputs of the function the names of the models to use to get the biomass compostion')
elseif sameBiomass == 0 && nargin == 4
    [~,idx_biomass] = ismember(id,names4biomass);
end

%% load data
folderNames = dir;
for i = 1:length(folderNames)
    names{i,1} = folderNames(i).name;
end
% get only relevant results
names = names(contains(names,'D0_Biosynthesis'));
%names = names(~contains(names,'pyruvate'));
% get results always on the same order
names = unique(names);
id = names;

for i = 1:length(names)
    myFolder = strcat('./',names{i},'/Models/',id{i},'/');
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    if ~isempty(matFiles)
        var = strrep(matFiles.name,'.mat','');
        info = load(strcat(matFiles.folder,'/',matFiles.name));
        info = info.(var);
        info = info.info_LMDP;
        %assignin('base',strcat('lump_info',num2str(i)),info);
        lump_info{i} = info;
    else
        lump_info{i} = [ ];
    end
end

for i = 1:length(names)
    myFolder = strcat('./',names{i},'/WorkSpaces/',id{i},'/',id{i},'/');
    filePattern = fullfile(myFolder, '*addBIOMASSnis_6.mat');
    matFiles = dir(filePattern);
    if ~isempty(matFiles)
        var = strrep(matFiles.name,'.mat','');
        moreInfo = load(strcat(matFiles.folder,'/',matFiles.name),'LumpCstModel','DPsAll','bbb_metnames','summaryOfLUMPS');
        moreInfo_lump{i} = moreInfo;
    else
        lump_info{i} = [ ];
    end
end

if ~idx_biomass
    % for carveMe models all models have the same biomass composition
    % so we get the info from one GEM
    i = 1;
    myFolder = strcat('./',names{i},'/Models/',strrep(id{i},'_D0_Biosynthesis',''),'/');
    filePattern = fullfile(myFolder, '*.mat');
    matFiles = dir(filePattern);
    var = strrep(matFiles.name,'.mat','');
    info = load(strcat(matFiles.folder,'/',matFiles.name));
    info = info.(var);
    commonBBBs = info.mets(find(info.S(:,find(info.c))<0));
else
    % in case we use models with different biomass compositions
    % we get all bbbs
    for k = 1:length(idx_biomass)
        myFolder = strcat('./',names{idx_biomass(k)},'/Models/',id{idx_biomass(k)},'/');
        filePattern = fullfile(myFolder, '*.mat');
        matFiles = dir(filePattern);
        var = strrep(matFiles.name,'.mat','');
        info = load(strcat(matFiles.folder,'/',matFiles.name));
        info = info.(var);
        info = info.mets(find(info.S(:,find(info.c))<0));
        assignin('base',strcat('commonBBBs',num2str(k)),info);
    end
    % and we keep the common ones
    commonBBBs = eval(strcat('commonBBBs',num2str(idx_biomass(1))));
    for k = 2:length(idx_biomass)
        bbbs = eval(strcat('commonBBBs',num2str(idx_biomass(k))));
        commonBBBs = commonBBBs(ismember(commonBBBs,bbbs));
    end
end

%%
% we remove some BBBs from the analysis
bbbsToExclude = {'atp_c';'h2o_c';'cl_c';'fe2_c';'so4_c';'k_c';'cu2_c';'ca2_c';'mg2_c';'mn2_c';'cobalt2_c';'zn2_c';'fe3_c'};
% to add the KEGG and SEED IDs
commonBBBs = setdiff(commonBBBs,bbbsToExclude);
% we can put BBBs into categories
% here we provide different categories and then we see which of the
% BBBs we have in the commonBBBs list
amino_acids = {'ala__L_c';'arg__L_c';'asn__L_c';'asp__L_c';'cys__L_c';'gln__L_c';'glu__L_c';'gly_c';'his__L_c';'ile__L_c';'leu__L_c';'lys__L_c';'met__L_c';'phe__L_c';'pro__L_c';'ser__L_c';'thr__L_c';'trp__L_c';'tyr__L_c';'val__L_c'};
amino_acids = [amino_acids;{'C00041_c';'C00062_c';'C00152_c';'C00049_c';'C00097_c';'C00064_c';'C00025_c';'C00037_c';'C00135_c';'C00407_c';'C00123_c';'C00047_c';'C00073_c';'C00079_c';'C00148_c';'C00065_c';'C00188_c';'C00078_c';'C00082_c';'C00183_c'}];
amino_acids = [amino_acids;{'cpd00035_c';'cpd00051_c';'cpd00132_c';'cpd00041_c';'cpd00084_c';'cpd00053_c';'cpd00023_c';'cpd00033_c';'cpd00119_c';'cpd00322_c';'cpd00107_c';'cpd00039_c';'cpd00060_c';'cpd00066_c';'cpd00129_c';'cpd00054_c';'cpd00161_c';'cpd00065_c';'cpd00069_c';'cpd00156_c'}];
amino_acids = amino_acids(ismember(amino_acids,commonBBBs));
nucleotides = {'ctp_c';'datp_c';'dctp_c';'dgtp_c';'dttp_c';'gtp_c';'utp_c'};
nucleotides = [nucleotides;{'C00063_c';'C00131_c';'C00458_c';'C00286_c';'C00459_c';'C00044_c';'C00075_c'}];
nucleotides = [nucleotides;{'cpd00052_c';'cpd00115_c';'cpd00356_c';'cpd00241_c';'cpd00357_c';'cpd00038_c';'cpd00062_c'}];
nucleotides = nucleotides(ismember(nucleotides,commonBBBs));
lipids = {'uaagmda_c';'C05898_c';'cpd03495_c';'cpd15531_c';'C06025_c';'C00344_c';'C00641_c';'C01595_c';'C01885_c';'C04046_c';'C05752_c';'C08362_c'};
lipids = lipids(ismember(lipids,commonBBBs));
cofactors = {'10fthf_c';'amet_c';'coa_c';'fad_c';'mlthf_c';'mql8_c';'nad_c';'nadp_c';'pydx5p_c';'ribflv_c';'thf_c';'thmpp_c'};
cofactors = [cofactors;{'C00828_c';'C00748_c';'C00032_c';'C00234_c';'C00019_c';'C00010_c';'C00016_c';'C00143_c';'C05819_c';'C00003_c';'C00006_c';'C00018_c';'C00255_c';'C00101_c';'C00068_c';'C00120_c';'C00440_c';'C00083_c';'C00129_c'}];
cofactors = [cofactors;{'cpd00201_c';'cpd00017_c';'cpd00010_c';'cpd00015_c';'cpd00125_c';'cpd15499_c';'cpd00003_c';'cpd00006_c';'cpd00016_c';'cpd00220_c';'cpd00087_c';'cpd00056_c'}];
cofactors = cofactors(ismember(cofactors,commonBBBs));
carbohydrates = {'glycogen_c';'C00182_c';'cpd00155_c';'cpd15514_c'};
carbohydrates = carbohydrates(ismember(carbohydrates,commonBBBs));

all_bbbs = {amino_acids,nucleotides,cofactors,lipids,carbohydrates};
all_bbbs_aux = {'amino_acids','nucleotides','cofactors','carbohydrates'};
categories = {'Amino_Acids';'Nucleotides';'Cofactors';'Lipids';'Carbohydrates'};
idx_empty = find(cellfun(@isempty,all_bbbs)==0);
categories = categories(idx_empty);
%% get statistics
all_c = {};
all_alts = [];
all_sizeAlt = [];
meanSize = [];
for k = 1:length(categories)
    U = {};
    sizeAlt = [];
    alts = [];
    commonBBBs1 = all_bbbs{k};

    for idx = 1:length(names)
        for idx2 = 1:length(commonBBBs1)
            info1 = lump_info{idx};
            if ~isempty(info1)
                f1 = find(ismember(info1(:,2),commonBBBs1{idx2}));
                s1 = [];
                for ss = 1:length(f1)
                    aux = info1(f1(ss,1),4);
                    s1(ss) = length(aux{1});
                end
                ms1 = min(s1);
                f1 = f1(ismember(s1,ms1));

                if not(isempty(f1))
                    alts(idx,idx2) = length(f1);
                    sizeAlt(idx,idx2) = length(info1{f1(1),4});
                    aux = info1(f1,4);
                    U{idx,idx2} = length(unique(vertcat(aux{:})));
                else
                    alts(idx,idx2) = 0;
                    sizeAlt(idx,idx2) = 0;
                    U{idx,idx2} = 0;

                end
            else
                alts(idx,idx2) = 0;
                sizeAlt(idx,idx2) = 0;
                U{idx,idx2} = 0;
            end

        end
    end
    all_c = [all_c  U];
    all_alts = [all_alts alts];
    all_sizeAlt = [all_sizeAlt sizeAlt];
    meanSize = [meanSize ;mean(sizeAlt,2)'];
end

% plot
yvalues = vertcat(all_bbbs{:});
%yvalues = vertcat(commonBBBs);
[~,idx] = ismember(yvalues,info.mets);
yvalues = info.metNames(idx); %to get the names instead of ids
figure(1)
h = heatmap(id,yvalues,all_alts','FontName','Avenir Next');
h.Title = 'Alternatives';
h.XLabel = 'Organisms';
h.YLabel = 'BBBs';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_blue;
caxis([0,max(all_alts(:))])
saveas(h,'Alternatives','svg')
T = array2table(all_alts', 'VariableNames', id, 'RowNames', yvalues);
filename = strcat('statistics_biosynthesis_Alternatives.xlsx');
writetable(T, filename, 'WriteRowNames', true);

figure(2)
h = heatmap(id,yvalues,all_sizeAlt','FontName','Avenir Next');
h.Title = 'Subnetwork size';
h.XLabel = 'Organisms';
h.YLabel = 'BBBs';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_red;
caxis([0,max(all_sizeAlt(:))])
saveas(h,'Subnetwork size','svg')
T = array2table(all_sizeAlt', 'VariableNames', id, 'RowNames', yvalues);
filename = strcat('statistics_biosynthesis_subnetwork_size.xlsx');
writetable(T, filename, 'WriteRowNames', true);


cl = clustergram(meanSize);
r = cell2dbl(cl.RowLabels);
c = cell2dbl(cl.ColumnLabels);

figure(3)
h = heatmap(id(c),categories(r),meanSize(r,c),'FontName','Avenir Next');
h.Title = 'Subnetwork mean size';
h.XLabel = 'Organisms';
h.YLabel = 'BBB category';
s = struct(h);
s.XAxis.TickLabelRotation = 90;
h.Colormap = map_red;
caxis([0,max(meanSize(:))])
saveas(h,'Subnetwork mean size','svg')

T = array2table(meanSize, 'VariableNames', id, 'RowNames', categories);
filename = strcat('statistics_biosynthesis_meanSize.xlsx');
writetable(T, filename, 'WriteRowNames', true);
%% comparison

for k = 1:length(categories)
    J_lumps_all ={};
    J_lumps = [];
    commonBBBs1 = all_bbbs{k};
    for idx = 1:length(names)
        for idx2 = 1:length(names)
            J_lumps = 0;
            info1 = lump_info{idx};
            info2 = lump_info{idx2};
            if ~isempty(info1) &&~isempty(info2)
                for i = 1 : length(commonBBBs1)

                    f1 = find(ismember(info1(:,2),commonBBBs1{i}));
                    f2 = find(ismember(info2(:,2),commonBBBs1{i}));

                    if not(isempty(f1)) && not(isempty(f2))
                        aux1 = info1(f1(:),4);
                        aux1 = unique(vertcat(aux1{:}));

                        aux2 = info2(f2(:),4);
                        aux2 = unique(vertcat(aux2{:}));

                        overlap = numel(find(ismember(aux1,aux2)));
                        union = numel(aux1) + numel(aux2) - overlap;
                        if jaccard
                            J_lumps(i,1) = overlap./union;
                        else
                            J_lumps(i,1) = overlap./numel(aux1) ;
                        end
                        clear overlap union aux1 aux2
                    else
                        J_lumps(i,1) = 0;
                    end
                end
            end
            J_lumps_all{idx,idx2} = J_lumps;
            J_average(idx,idx2) = mean(J_lumps);

        end
        J_BBBs{k} = J_lumps_all;
        J_BBBs_average{k} = J_average;

    end
end
%save all data to a D0_Biosynthesis file
save('comparison_bioynthesis.mat','J_BBBs','J_BBBs_average', 'id','categories')

%save average in xlsx for visualization and already plot a heatmap
for i = 1:length(categories)
    aux = J_BBBs_average{i};
    if ~isempty(aux)


        cl = clustergram(aux);
        r = cell2dbl(cl.RowLabels);

        figure(3+i)
        h = heatmap(id(r),id(r),aux(r,r),'FontName','Avenir Next');
        h.Title = categories{i};
        h.XLabel = 'Organisms';
        h.YLabel = 'Organisms';
        s = struct(h);
        s.XAxis.TickLabelRotation = 90;
        h.Colormap = map_red;
        saveas(h,strcat('average_similairity_',categories{i}),'svg')

        T = array2table(aux, 'VariableNames', id, 'RowNames', id);
        filename = strcat('average_similarity_new_',categories{i},'.xlsx');
        writetable(T, filename, 'WriteRowNames', true);
    end
end


%% extract the cost of biosynthesis for each one of the biomass building blocks for each organism
for idx = 1: length(names)
     model_all = createModel();
    cost = {};
    info1 = lump_info{idx};
    if ~isempty(info1)
    model1_all = createModel();
    for i = 1 : length(commonBBBs)

        f1 = find(ismember(info1(:,2),commonBBBs{i}));

        if not(isempty(f1))

            rxns1 = info1(f1,1);
            formulas1 = info1(f1,3);
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
            coeff = [];
%             aux = strsplit(commonBBBs{i},'_');
%             aux = aux{1};
            f = find(ismember(model1.mets,commonBBBs{i}));
            if f
                coeff = full(model1.S(f,:));
                model1.S = model1.S./coeff;
            end
            if not(ismember(commonBBBs{i},model1.mets)) && isfield(more,'DPsAll')
                more1 = moreInfo_lump{idx};
                f = find(ismember(more1.LumpCstModel.mets,commonBBBs{i}));
                idxNo = [find(more1.LumpCstModel.isTrans);  find(contains(more1.LumpCstModel.rxns,'EX_'))]; %need to fix this
                % get drains and transporters
                S = find(more1.LumpCstModel.S(f,:));
                S = setdiff(S,idxNo);
                f = find(ismember(more1.bbb_metnames,commonBBBs{i}));
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
                    model1.rxns{end+1} = strcat('aux_',commonBBBs{i}); model1.lb(end+1) = -100; model1.ub(end+1) = 100;
                    if model1.S(find(ismember(model1.mets,commonBBBs{i})),end)<0
                        model1.S(:,end) = -1*model1.S(:,end);
                    end
                    model1 = removeRxns(model1,rxns);
                end
            end
            coeff1 = sum(model1.S,2);
            coeff1 = coeff1/length(rxns1);
            model1.S(:,end+1) = coeff1;
            model1.rxns{end+1} = commonBBBs{i};
            cost(i,1) = printRxnFormulaCost(model1,model1.rxns(end),0);
            model1_all = addReactionTFA(model1_all,commonBBBs{i},cost{end});


        else
            cost(i,1) = {':'};
        end
    end
    T = array2table(cost, 'RowNames', commonBBBs);
    filename = strcat('cost_new_',id{idx},'.xlsx');
    writetable(T, filename, 'WriteRowNames', true);
    

    all_costs{idx,1} = cost;
%     T = array2table(cost, 'RowNames', Tasks);
%     filename = strcat('cost_',id{idx},'.xlsx');
%     writetable(T, filename, 'WriteRowNames', true);
    for t = 1:length(commonBBBs)
        model_all = addReactionTFAnis(model_all,strcat(id{idx},'_',commonBBBs{t}),cost{t,1});
    end
    all_rxns{idx} = model_all;
end
end
%% check if we have the same BBBs in the average lump
for t = 1:length(commonBBBs)
    J_cost = [];
    for idx = 1:length(id)
        
        model1 = all_rxns{idx};
        f1 = find(ismember(model1.rxns,strcat(id{idx},'_',commonBBBs(t))));
        m1 = find(model1.S(:,f1));

        for idx2 = 1:length(id)
           
            model2 = all_rxns{idx2};
            f2 = find(ismember(model2.rxns,strcat(id{idx2},'_',commonBBBs(t))));

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

save('comparison_bioynthesis_cost.mat','J_Tasks_cost', 'id','commonBBBs')
id = strrep(id,'_','-');
Tasks = strrep(commonBBBs,'_','-');
%save average in xlsx for visualization and already plot a heatmap
for i = 1:length(Tasks)
    if i~=54 && i~=55 && i~=15&& i~=23&&i~=35&&i~=42&&i~=46
    aux = J_Tasks_cost{i,1};
    if ~isempty(aux)

aux(isnan(aux)) = 0;

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
end
%%  now just get the coeff for ATP, NADH, NADPH, NAD and NADP
aux = [];
mets = {'PHB_c';'cpd00002_c';'cpd00008_c';'C00003_c';'C00004_c';'C00006_c';'C00005_c'};
metNames = {'PHB';'ATP';'ADP';'NAD';'NADH' ;'NADP';'NADPH'};
Tasks = strrep(Tasks,'-','_');
id = strrep(id,'-','_');
for t = 1:length(Tasks)
    coeff = [];
    for idx = 1:length(id)
       % if idx~=7 && idx~=15
        model1 = all_rxns{idx};
        f1 = find(ismember(model1.rxns,strcat(id(idx),'_',Tasks(t))));
        [~,fm] = ismember(mets,model1.mets);
        for i = 1:length(mets)
            if fm(i)
                aux(i,1) = full(model1.S(fm(i),f1));
            else
                aux(i,1) = 0;
            end
        end
        coeff= [coeff ; aux'];
%         else
%           coeff= [coeff ; zeros(length(mets),1)'];  
%   
%         end
    end
          Coeff_Tasks{t,1} = coeff;
    
end

save('comparison_bioynthesis_costCoeff.mat','Coeff_Tasks', 'mets','id','metNames')
id = strrep(id,'_','-');
Tasks = strrep(Tasks,'_','-');
%save average in xlsx for visualization and already plot a heatmap
clear aux
for i = 1:length(Tasks)
    aux = Coeff_Tasks{i,1};
    if ~isempty(aux)

%         cl = clustergram(aux);
%         r = cell2dbl(cl.RowLabels);

        figure()
        aux = -aux;
        for k = 1:size(aux,2)

            f = find(aux(:,k) == -0);
            aux(f,k) = 0;
        end
        %h = heatmap(metNames,id(r),aux(r,:),'FontName','Avenir Next');
        h = heatmap(metNames,id,aux,'FontName','Avenir Next');
        h.Title = Tasks{i};
        %h.XLabel = 'Organisms';
        h.YLabel = 'Organisms';
        s = struct(h);
        s.XAxis.TickLabelRotation = 90;
        h.Colormap = map_red;
        saveas(h,strcat('costCoeff_',Tasks{i}),'svg')

        T = array2table(-aux, 'VariableNames', metNames, 'RowNames', id);
        filename = strcat('costCoeff_',Tasks{i},'.xlsx');
        writetable(T, filename, 'WriteRowNames', true);
    end
end