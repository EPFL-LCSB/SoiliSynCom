% This function modifies the metabolic model to fix faulty substrate cycles.
% The function finds reactions involving substrate ('C00007_c') as a reactant
% or product, and modifies the lower and upper bounds of those reactions.
% The modified reactions, their formulas and fluxes are then written to an
% excel file. The function returns the final task report, the modified model,
% and two lists of reactions involving substrate as a reactant ('rxnsub') or
% product ('rxnslb').
function modifiedModel = FixFalsePositives(model,excelFileName)

    % Set the optimization solver to 'gurobi'
    setRavenSolver('gurobi');
    if (nargin<4)
        excelFileName = 'ATP_ResponsibleRxns.xlsx';
    end
    
    % Run the task 'Tasks_substrate.xlsx' to get the essential reactions and fluxes
    [verificationTaskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,'Tasks.xlsx',true,true,true);
    
    % Define the ID for substrate
    substrate = 'C00002';
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
        
        % Get the active reactions and their fluxes
        idx = find(essentialFluxes);
        ActiveRxns = model.rxns(idx);
        fluxes=essentialFluxes(idx);
        
        % Define the reactions to be excluded
        excludedReactions = {'EXC_BOTH', 'T_p_to_c', 'transport','T_e_to_p','Biomass','T_c_to_e','ATPM','Rnew','R09077'};
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
        
        
        % Iterate through each reaction involving substrate
        for i = 1:length(Formulas)
            
            % Split the reaction into reactants and products
            split_reaction = regexp(Formulas{i}, ' <=> | -> | ,-', 'split');
            reactants = split_reaction{1};
            products = split_reaction{2};
            
            % If substrate is a reactant, add the reaction to 'rxnslb'
            if fluxes{i}<-20 && fluxes{i}>-50
                rxnslb{end + 1, 1} = ActiveRxns{i};
                formulasSlb{end + 1, 1} = Formulas{i};
                fluxeslb(end + 1, 1) = fluxes(i);
            end
            
            % If substrate is a product, add the reaction to 'rxnsub'
            if fluxes{i}>20 && fluxes{i}<50
                rxnsub{end + 1, 1} = ActiveRxns{i};
                formulasSub{end + 1, 1} = Formulas{i};
                fluxesub(end + 1, 1) = fluxes(i);
            end
        end
        % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
        f = find(ismember(model.rxns, rxnsub));
        model.ub(f) = 0;
        
        % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
        f = find(ismember(model.rxns, rxnslb));
        model.lb(f) = 0;
        % Run the task 'Tasks_substrate.xlsx' again to get the new essential reactions and fluxes
        [verificationTaskReport, essentialRxns, taskStructure, essentialFluxes]=checkTasks(model,'Tasks.xlsx',true,true,true);
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
    writetable(T, excelFileName, 'WriteRowNames', true);
    
    % Display a message indicating that the Excel file has been created
    display("The reactions modified with their details are saved in the excel file 'substrate containing reactions modified.xlsx'")
end
