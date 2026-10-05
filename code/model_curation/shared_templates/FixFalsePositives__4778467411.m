% This function modifies the metabolic model to fix faulty substrate cycles.
% The function finds reactions involving substrate ('C00007_c') as a reactant
% or product, and modifies the lower and upper bounds of those reactions.
% The modified reactions, their formulas and fluxes are then written to an
% excel file. The function returns the final task report, the modified model,
% and two lists of reactions involving substrate as a reactant ('rxnsub') or
% product ('rxnslb').
function modifiedModel = FixFalsePositives(model,verificationTasks,Substrates,excelFileName)

    % Set the optimization solver to 'gurobi'
    setRavenSolver('gurobi');
    if (nargin<4)
        excelFileName = 'FP_ResponsibleRxns.xlsx';
    end
    SpecificSubstrates={'C00042';'C00109';'C00093';'C00103';'C00140';'C03383';'C00065';'C00134';'C00156';'C00334'};

    for v=1:numel(Substrates)
        % Run the task 'Tasks_substrate.xlsx' to get the essential reactions and fluxes
        [verificationTaskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,[],true,true,true,verificationTasks(v));
        
        % Define the ID for substrate
        substrate = Substrates(v);
        substrateName=model.metNames(find(ismember(model.mets,strcat(substrate,'_c'))));
        % Initialize lists for storing reaction IDs, formulas and fluxes
        rxnsub = {}; % Reactions with substrate as a product
        rxnslb = {}; % Reactions with substrate as a reactant
        formulasSub = {};
        formulasSlb = {};
        fluxesub = {};
        fluxeslb = {};
        
        % Repeat the following steps until no tasks are OK (taskReport.ok = 0)
        while verificationTaskReport.ok ~= 0
            if any(ismember(SpecificSubstrates,substrate))
                
                % Get the active reactions and their fluxes
                idx = find(essentialFluxes);
                ActiveRxns = model.rxns(idx);
                fluxes=essentialFluxes(idx);
                
                % Define the reactions to be excluded
                excludedReactions = {'EXC_BOTH', 'T_p_to_c', 'transport','T_e_to_p','Biomass','T_c_to_e'};
                subSystems = model.subSystems(find(ismember(model.rxns, ActiveRxns)));
                % Loop over the ActiveRxns and check if they contain any of the excluded reactions
                ActiveRxnsFiltered = {};
                fluxesFiltered={};
                for i = 1:length(ActiveRxns)
                    if ~any(contains(ActiveRxns(i), excludedReactions)) && any(~strcmp(subSystems{i}, 'GapFilled'))
                        ActiveRxnsFiltered(end+1) = ActiveRxns(i);
                        fluxesFiltered(end+1)=num2cell(fluxes(i));
                    end
                end
                % Get the formulas of the active reactions
                ActiveRxns = ActiveRxnsFiltered;
                fluxes=fluxesFiltered;
                Formulas = printRxnFormula(model, ActiveRxns);
                
                % Find the reactions that involve substrate
                matches = regexp(Formulas, substrate, "match");
                idxsubstrate = find(~cellfun('isempty', matches));
                
                % Get the active reactions involving substrate, their formulas, and their fluxes
                ActiveRxnssubstrate = ActiveRxns(idxsubstrate);
                Formulassubstrate = printRxnFormula(model, ActiveRxnssubstrate);
                fluxessubstrate = fluxes(idxsubstrate);
                
                % Iterate through each reaction involving substrate
                for i = 1:length(Formulassubstrate)
                    
                    % Split the reaction into reactants and products
                    split_reaction = regexp(Formulassubstrate{i}, ' <=> | -> | ,-', 'split');
                    reactants = split_reaction{1};
                    products = split_reaction{2};
                    
                    % If substrate is a reactant, add the reaction to 'rxnslb'
                    if contains(reactants, substrate) && fluxessubstrate{i}>20
                        rxnslb{end + 1, 1} = ActiveRxnssubstrate{i};
                        formulasSlb{end + 1, 1} = Formulassubstrate{i};
                        fluxeslb(end + 1, 1) = fluxessubstrate(i);
                    end
                    
                    % If substrate is a product, add the reaction to 'rxnsub'
                    if contains(products, substrate) && fluxessubstrate{i}<-20
                        rxnsub{end + 1, 1} = ActiveRxnssubstrate{i};
                        formulasSub{end + 1, 1} = Formulassubstrate{i};
                        fluxesub(end + 1, 1) = fluxessubstrate(i);
                    end
                end
                
                % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
                f = find(ismember(model.rxns, rxnsub));
                model.lb(f) = 0;
                
                % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
                f = find(ismember(model.rxns, rxnslb));
                model.ub(f) = 0;        

                % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
                f = find(ismember(model.rxns, rxnsub));
                model.lb(f) = 0;
                
                % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
                f = find(ismember(model.rxns, rxnslb));
                model.ub(f) = 0;
                
                % Run the task 'Tasks_substrate.xlsx' again to get the new essential reactions and fluxes
                [verificationTaskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,[],true,true,true,verificationTasks(v));
            else
                
                % Get the active reactions and their fluxes
                idx = find(essentialFluxes);
                ActiveRxns = model.rxns(idx);
                fluxes=essentialFluxes(idx);
                
                % Define the reactions to be excluded
                excludedReactions = {'EXC_BOTH', 'T_p_to_c', 'transport','T_e_to_p','Biomass','T_c_to_e'};
                subSystems = model.subSystems(find(ismember(model.rxns, ActiveRxns)));
                % Loop over the ActiveRxns and check if they contain any of the excluded reactions
                ActiveRxnsFiltered = {};
                fluxesFiltered={};
                for i = 1:length(ActiveRxns)
                    if ~any(contains(ActiveRxns(i), excludedReactions)) && any(~strcmp(subSystems{i}, 'GapFilled'))
                        ActiveRxnsFiltered(end+1) = ActiveRxns(i);
                        fluxesFiltered(end+1)=num2cell(fluxes(i));
                    end
                end
                % Get the formulas of the active reactions
                ActiveRxns = ActiveRxnsFiltered;
                fluxes=fluxesFiltered;
                Formulas = printRxnFormula(model, ActiveRxns);
                
                % Find the reactions that involve substrate
                matches = regexp(Formulas, substrate, "match");
                idxsubstrate = find(~cellfun('isempty', matches));
                
                % Get the active reactions involving substrate, their formulas, and their fluxes
                ActiveRxnssubstrate = ActiveRxns(idxsubstrate);
                Formulassubstrate = printRxnFormula(model, ActiveRxnssubstrate);
                fluxessubstrate = fluxes(idxsubstrate);
                
                % Iterate through each reaction involving substrate
                for i = 1:length(Formulassubstrate)
                    
                    % Split the reaction into reactants and products
                    split_reaction = regexp(Formulassubstrate{i}, ' <=> | -> | ,-', 'split');
                    reactants = split_reaction{1};
                    products = split_reaction{2};
                    
                    % If substrate is a reactant, add the reaction to 'rxnslb'
                    if contains(reactants, substrate)
                        rxnslb{end + 1, 1} = ActiveRxnssubstrate{i};
                        formulasSlb{end + 1, 1} = Formulassubstrate{i};
                        fluxeslb(end + 1, 1) = fluxessubstrate(i);
                    end
                    
                    % If substrate is a product, add the reaction to 'rxnsub'
                    if contains(products, substrate)
                        rxnsub{end + 1, 1} = ActiveRxnssubstrate{i};
                        formulasSub{end + 1, 1} = Formulassubstrate{i};
                        fluxesub(end + 1, 1) = fluxessubstrate(i);
                    end
                end
                
                % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
                f = find(ismember(model.rxns, rxnsub));
                model.lb(f) = 0;
                
                % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
                f = find(ismember(model.rxns, rxnslb));
                model.ub(f) = 0;
                
                % Run the task 'Tasks_substrate.xlsx' again to get the new essential reactions and fluxes
                [verificationTaskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,[],true,true,true,verificationTasks(v));
            end
        end
        % Assign the final modified model to 'modifiedModel'
        modifiedModel = model;
        
        % Merge the lists of reaction IDs, formulas, and fluxes
        rxnall = [rxnsub; rxnslb];
        formulasAll = [formulasSub; formulasSlb];
        fluxAll = [fluxesub; fluxeslb];
        rxnNames=printRxnFormula(model,rxnall,1,1,1);
        % Display the total number of modified reactions
        fprintf('Number of substrate containing reactions modified: %d\n', length(rxnall))
        
        % Write the modified reactions, their formulas, and their fluxes to an Excel file
        rxnFormulasFluxes = [rxnall, formulasAll, rxnNames,num2cell(fluxAll)];
        T = cell2table(rxnFormulasFluxes, 'VariableNames', {'Reaction', 'Formula','Names','Flux'});
        substrateName=erase(substrateName,'-');
        substrateName=erase(substrateName,'/');
        writetable(T, excelFileName, 'Sheet', cell2str(strcat(substrate,'-',substrateName)), 'WriteRowNames', true);
        
        % Display a message indicating that the Excel file has been created
        display("The reactions modified with their details are saved in the excel file 'substrate containing reactions modified.xlsx'")
    end
end
