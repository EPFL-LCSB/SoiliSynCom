% This function modifies the metabolic model to fix faulty oxygen cycles.
% The function finds reactions involving oxygen ('C00007_c') as a reactant
% or product, and modifies the lower and upper bounds of those reactions.
% The modified reactions, their formulas and fluxes are then written to an
% excel file. The function returns the final task report, the modified model,
% and two lists of reactions involving oxygen as a reactant ('rxnsub') or
% product ('rxnslb').
function modifiedModel = FixFaultyOxygenCycles(model)

    % Initialize reaction bounds for 'EXC_BOTH_C00007_e' in the model
    f = find(contains(model.rxns,'EXC_BOTH_C00007_e'));
    model.lb(f) = 0;
    model.ub(f) = 0;

    % Set the optimization solver to 'gurobi'
    setRavenSolver('gurobi');

    % Run the task 'Tasks_Oxygen.xlsx' to get the essential reactions and fluxes
    [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks_Oxygen.xlsx', true, 'getEssential', true);
     % Check if taskReport is false, if it is then print a message and end the function
    if taskReport.ok == 0
        fprintf('No faulty oxygen cycles found.\n');
            f = find(contains(model.rxns,'EXC_BOTH_C00007_e'));
            model.lb(f) = -20;
            model.ub(f) = 50;
        modifiedModel = model;
        rxnsub = {};
        rxnslb = {};
        return
    end
    % Define the ID for oxygen
    oxygen = 'C00007_c';

    % Initialize lists for storing reaction IDs, formulas and fluxes
    rxnsub = {}; % Reactions with oxygen as a product
    rxnslb = {}; % Reactions with oxygen as a reactant
    formulasSub = {};
    formulasSlb = {};
    fluxesub = [];
    fluxeslb = [];

    % Repeat the following steps until no tasks are OK (taskReport.ok = 0)
    while taskReport.ok ~= 0

        % Get the active reactions and their fluxes
        idx = find(essentialFluxes);
        ActiveRxns = model.rxns(idx);
        fluxes=essentialFluxes(idx);
        if ismember('R09077', ActiveRxns)
            idx_to_exclude = find(ActiveRxns == "R09077");
            ActiveRxns(idx_to_exclude) = [];
            fluxes(idx_to_exclude) = [];
        end
        
        % Get the formulas of the active reactions
        Formulas = printRxnFormula(model, ActiveRxns);

        % Find the reactions that involve oxygen
        matches = regexp(Formulas, oxygen, "match");
        idxOxy = find(~cellfun('isempty', matches));
        
        % Get the active reactions involving oxygen, their formulas, and their fluxes
        ActiveRxnsOxy = ActiveRxns(idxOxy);
        FormulasOxy = printRxnFormula(model, ActiveRxnsOxy);
        fluxesOxy = fluxes(idxOxy);

        % Iterate through each reaction involving oxygen
        for i = 1:length(FormulasOxy)

            % Split the reaction into reactants and products
            split_reaction = regexp(FormulasOxy{i}, ' <=> | -> | ,-', 'split');
            reactants = split_reaction{1};
            products = split_reaction{2};

            % If oxygen is a reactant, add the reaction to 'rxnslb'
            if contains(reactants, oxygen)
                rxnslb{end + 1, 1} = ActiveRxnsOxy{i};
                formulasSlb{end + 1, 1} = FormulasOxy{i};
                fluxeslb(end + 1, 1) = fluxesOxy(i);
            end
            
            % If oxygen is a product, add the reaction to 'rxnsub'
            if contains(products, oxygen)
                rxnsub{end + 1, 1} = ActiveRxnsOxy{i};
                formulasSub{end + 1, 1} = FormulasOxy{i};
                fluxesub(end + 1, 1) = fluxesOxy(i);
            end
        end

        % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
        f = find(ismember(model.rxns, rxnsub));
        model.ub(f) = 0;
        
        % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
        f = find(ismember(model.rxns, rxnslb));
        model.lb(f) = 0;

        % Run the task 'Tasks_Oxygen.xlsx' again to get the new essential reactions and fluxes
        [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks_Oxygen.xlsx', true, 'getEssential', true);
    end

    % Assign the final modified model to 'modifiedModel'
    f = find(contains(model.rxns,'EXC_BOTH_C00007_e'));
    model.lb(f) = -20;
    model.ub(f) = 50;
    
    modifiedModel = model;
    
    % Restore the bounds of 'EXC_BOTH_C00007_e' in the model
    f = find(contains(model.rxns,'EXC_BOTH_C00007_e'));
    model.lb(f) = -20;
    model.ub(f) = 50;

    % Merge the lists of reaction IDs, formulas, and fluxes
    rxnall = [rxnsub; rxnslb];
    formulasAll = [formulasSub; formulasSlb];
    fluxAll = [fluxesub; fluxeslb];
    
    % Display the total number of modified reactions
    fprintf('Number of Oxygen containing reactions modified: %d\n', length(rxnall))
    
    % Write the modified reactions, their formulas, and their fluxes to an Excel file
    rxnFormulasFluxes = [rxnall, formulasAll, num2cell(fluxAll)];
    T = cell2table(rxnFormulasFluxes, 'VariableNames', {'Reaction', 'Formula', 'Flux'});
    writetable(T, 'Oxygen containing reactions modified.xlsx');
    % Display a message indicating that the Excel file has been created
    display("The reactions modified with their details are saved in the excel file 'Oxygen containing reactions modified.xlsx'")
end
