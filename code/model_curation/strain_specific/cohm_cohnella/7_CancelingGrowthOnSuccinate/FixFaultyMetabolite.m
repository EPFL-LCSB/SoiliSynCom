% This function modifies the metabolic model to fix faulty FaultyMetabolite cycles.
% The function finds reactions involving FaultyMetabolite ('C00007_c') as a reactant
% or product, and modifies the lower and upper bounds of those reactions.
% The modified reactions, their formulas and fluxes are then written to an
% excel file. The function returns the final task report, the modified model,
% and two lists of reactions involving FaultyMetabolite as a reactant ('rxnsub') or
% product ('rxnslb').
function modifiedModel = FixFaultyMetabolite(model,Metabolite)
    % Set the optimization solver to 'gurobi'
    setRavenSolver('gurobi');

    % Run the task 'Tasks_FaultyMetabolite.xlsx' to get the essential reactions and fluxes
    [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks.xlsx', true, 'getEssential', true);
    
    % Define the ID for FaultyMetabolite
    FaultyMetabolite = Metabolite;

    % Initialize lists for storing reaction IDs, formulas and fluxes
    rxnsub = {}; % Reactions with FaultyMetabolite as a product
    rxnslb = {}; % Reactions with FaultyMetabolite as a reactant
    formulasSub = {};
    formulasSlb = {};
    fluxesub = {};
    fluxeslb = {};

    % Repeat the following steps until no tasks are OK (taskReport.ok = 0)
    while taskReport.ok ~= 0

        % Get the active reactions and their fluxes
        idx = find(essentialFluxes);
        ActiveRxns = model.rxns(idx);
        fluxes=essentialFluxes(idx);
        
        % Define the reactions to be excluded
        excludedReactions = {'EXC_BOTH', 'T_p_to_c', 'transport'};
        
        % Loop over the ActiveRxns and check if they contain any of the excluded reactions
        ActiveRxnsFiltered = {};
        fluxesFiltered={};
        for i = 1:length(ActiveRxns)
            if ~any(contains(ActiveRxns(i), excludedReactions))
                ActiveRxnsFiltered(end+1) = ActiveRxns(i);
                fluxesFiltered(end+1)=num2cell(fluxes(i));
            end
        end
        % Get the formulas of the active reactions
        ActiveRxns = ActiveRxnsFiltered;
        fluxes=fluxesFiltered;
        Formulas = printRxnFormula(model, ActiveRxns);

        % Find the reactions that involve FaultyMetabolite
        matches = regexp(Formulas, FaultyMetabolite, "match");
        idxFaultyMetabolite = find(~cellfun('isempty', matches));
        
        % Get the active reactions involving FaultyMetabolite, their formulas, and their fluxes
        ActiveRxnsFaultyMetabolite = ActiveRxns(idxFaultyMetabolite);
        FormulasFaultyMetabolite = printRxnFormula(model, ActiveRxnsFaultyMetabolite);
        fluxesFaultyMetabolite = fluxes(idxFaultyMetabolite);

        % Iterate through each reaction involving FaultyMetabolite
        for i = 1:length(FormulasFaultyMetabolite)

            % Split the reaction into reactants and products
            split_reaction = regexp(FormulasFaultyMetabolite{i}, ' <=> | -> | ,-', 'split');
            reactants = split_reaction{1};
            products = split_reaction{2};

            % If FaultyMetabolite is a reactant, add the reaction to 'rxnslb'
            if contains(reactants, FaultyMetabolite)
                rxnslb{end + 1, 1} = ActiveRxnsFaultyMetabolite{i};
                formulasSlb{end + 1, 1} = FormulasFaultyMetabolite{i};
                fluxeslb(end + 1, 1) = fluxesFaultyMetabolite(i);
            end
            
            % If FaultyMetabolite is a product, add the reaction to 'rxnsub'
            if contains(products, FaultyMetabolite)
                rxnsub{end + 1, 1} = ActiveRxnsFaultyMetabolite{i};
                formulasSub{end + 1, 1} = FormulasFaultyMetabolite{i};
                fluxesub(end + 1, 1) = fluxesFaultyMetabolite(i);
            end
        end

        % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
        f = find(ismember(model.rxns, rxnsub));
        model.lb(f) = 0;
        
        % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
        f = find(ismember(model.rxns, rxnslb));
        model.ub(f) = 0;

        % Run the task 'Tasks_FaultyMetabolite.xlsx' again to get the new essential reactions and fluxes
        [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks.xlsx', true, 'getEssential', true);
    end

    % Assign the final modified model to 'modifiedModel'
    modifiedModel = model;
    
    % Merge the lists of reaction IDs, formulas, and fluxes
    rxnall = [rxnsub; rxnslb];
    formulasAll = [formulasSub; formulasSlb];
    fluxAll = [fluxesub; fluxeslb];
    
    % Display the total number of modified reactions
    fprintf('Number of FaultyMetabolite containing reactions modified: %d\n', length(rxnall))
    
    % Write the modified reactions, their formulas, and their fluxes to an Excel file
    rxnFormulasFluxes = [rxnall, formulasAll, num2cell(fluxAll)];
    T = cell2table(rxnFormulasFluxes, 'VariableNames', {'Reaction', 'Formula', 'Flux'});
    writetable(T, strcat(FaultyMetabolite,'containing reactions modified.xlsx'));
end
