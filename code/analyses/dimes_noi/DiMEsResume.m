addpath(genpath('./remind'));
addpath(genpath('${MATTFA_ROOT}'))%adding matTFA to matlab, i already
changeCobraSolver('cplex_direct');
saving_directory = './DiMEs_output/';
Models = {'ttmodel_bbet', 'ttmodel_buol', 'ttmodel_cauj', 'ttmodel_cohm', ...
    'ttmodel_cpi', 'ttmodel_cuvk', 'ttmodel_cxy', 'ttmodel_dri', 'ttmodel_flk', ...
    'ttmodel_lgu', 'ttmodel_lur', 'ttmodel_mgk', 'ttmodel_mlot', ...
    'ttmodel_pmor', 'ttmodel_pscf', 'ttmodel_pzu', 'ttmodel_racd', ...
    'ttmodel_rfa', 'ttmodel_tro', 'ttmodel_vpe', 'ttmodel_mip'};
targetNumAlt = 300;
carbonSource = {''};
% Check how many alternatives have already been computed
p= [1,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1];
p= 1./p; % percentages for parsimonious, here you can discretise as you wish
for m=12:17
model = eval(Models{m});
modeldescription = model.id;
model.description = modeldescription;
filename = strcat(modeldescription,'_PhenoMappingExchanges');

%%
    %parallel (i.e. different matlab instances)
    for i = 1:length(p)
        saving_directory_new = strcat(strcat(saving_directory,num2str(1/p(i)*100),'%_opt_yield_',cell2str(carbonSource),'/'));
        filename = strcat(modeldescription,'_PhenoMappingExchanges');

        r1 = load(strcat(saving_directory_new,filename,'.mat'));
        currentNumAlt = size(r1.DPsimm, 2);%here we can run sequentially for the different yield cuts but also in
        if currentNumAlt<150
           fprintf('Not Resuming DiMEs analysis: Found %d solutions\n', currentNumAlt);
            continue
        end        
        %printConstraint(model_p,{'CUT_pUP'})
        fprintf('Generating alternatives for %4.2f % of the optimal biomass yield \n',1/p(i)*100)
        
        sol = solveTFAmodelCplex(r1.model_p,r1.time);
        if ~isempty(r1.sol.val)
            
            [DPsimm_new,modelmilp2_new] = findDPMax_DiMEs(r1.modelmilp2, r1.NumAlt, r1.modelmilp2.indUSE,r1.time, r1.tagMin, strcat(saving_directory_new,filename),r1.tagMinSize);
            modelmilp2=modelmilp2_new;
            DPsimm=[r1.DPsimm DPsimm_new];
            % save workspace
            save(strcat(saving_directory_new,filename,'.mat'));
            
            % save results
            filename = strcat(modeldescription,'_PhenoMappingExchanges2');
            %r1 = load(strcat(saving_directory_new,filename,'.mat'));
            
            if not(isempty(DPsimm))
                % Extract info about composition of the IMMs
                [immOutput.StatsMets, immOutput.Mets, immOutput.drainClass, ...
                    immOutput.statsClass,immOutput.result,immOutput.Lengths,immOutput.Alt,Active] = extractInfo_DiMEs(r1.model_p, r1.DPsimm, ...
                    r1.model_p.indUSE,  r1.model_p.rowDrain_F,...
                    r1.model_p.rowDrain_R,r1.model_p.F_Flux,  r1.model_p.R_Flux) ;
                
                
                immOutput.summary = [strrep(immOutput.Mets(:,1),',',''), ...
                    num2cell(immOutput.result),num2cell(immOutput.Lengths),num2cell(immOutput.Alt),num2cell(immOutput.StatsMets)];
                
                writeData(strcat(saving_directory_new,filename,'_',num2str(1/p(i)*100),'.csv'), immOutput.summary,...
                    '%s\t%i\t%i\t%i\t%i', {'Drain ID', ...
                    'P/S/B','Lengths','Alt','Apperance in IMM'}, ...
                    '%s\t%s\t%s\t%s\t%s');
                
                save(strcat(saving_directory_new,'Active2_',modeldescription,'.mat'),'Active');
                clear Active
                
                elimIntermPMFile(saving_directory_new,filename)
            end
            f = find(contains(r1.model_p.constraintNames,'CUTimm'));
            r1.model_p.constraintNames(f) = strcat(r1.model_p.constraintNames(f),'_old');
        end
    end
    
    % Get unique solutions and prepare the input for the ILP step of ReMIND
    getUniqueSol(p,modeldescription,saving_directory)
    run GenerateSummaryTables.m
end