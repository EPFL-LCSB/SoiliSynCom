% This function modifies the metabolic model to fix faulty carbon cycles.
% The function finds reactions involving carbon ('C00007_c') as a reactant
% or product, and modifies the lower and upper bounds of those reactions.
% The modified reactions, their formulas and fluxes are then written to an
% excel file. The function returns the final task report, the modified model,
% and two lists of reactions involving carbon as a reactant ('rxnsub') or
% product ('rxnslb').
function modifiedModel = FixFaultyCarbonFixation(model)
    load('C02FixationTest.mat')
    load('PTYGData.mat')
    % Initialize reaction bounds for 'EXC_BOTH_C00007_e' in the model
    f=find(contains(model.rxns,EXC));
    model.lb(f)=0;
    model.ub(f)=0;
    f=find(contains(model.rxns,CarbonFixationTest));
    model.lb(f)=-20;
    model.ub(f)=50;

    % Set the optimization solver to 'gurobi'
    setRavenSolver('gurobi');

    % Run the task 'Tasks_carbon.xlsx' to get the essential reactions and fluxes
    [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks_Co2.xlsx', true, 'getEssential', true);
    
    % Define the ID for carbon
    carbon = 'C00011_c';

    % Initialize lists for storing reaction IDs, formulas and fluxes
    rxnsub = {}; % Reactions with carbon as a product
    rxnslb = {}; % Reactions with carbon as a reactant
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

        % Find the reactions that involve carbon
        matches = regexp(Formulas, carbon, "match");
        idxcarbon = find(~cellfun('isempty', matches));
        
        % Get the active reactions involving carbon, their formulas, and their fluxes
        ActiveRxnscarbon = ActiveRxns(idxcarbon);
        Formulascarbon = printRxnFormula(model, ActiveRxnscarbon);
        fluxescarbon = fluxes(idxcarbon);

        % Iterate through each reaction involving carbon
        for i = 1:length(Formulascarbon)

            % Split the reaction into reactants and products
            split_reaction = regexp(Formulascarbon{i}, ' <=> | -> | ,-', 'split');
            reactants = split_reaction{1};
            products = split_reaction{2};

            % If carbon is a reactant, add the reaction to 'rxnslb'
            if contains(reactants, carbon)
                rxnslb{end + 1, 1} = ActiveRxnscarbon{i};
                formulasSlb{end + 1, 1} = Formulascarbon{i};
                fluxeslb(end + 1, 1) = fluxescarbon(i);
            end
            
            % If carbon is a product, add the reaction to 'rxnsub'
            if contains(products, carbon)
                rxnsub{end + 1, 1} = ActiveRxnscarbon{i};
                formulasSub{end + 1, 1} = Formulascarbon{i};
                fluxesub(end + 1, 1) = fluxescarbon(i);
            end
        end

        % Modify the model by setting the upper bound of 'rxnsub' reactions to 0
        f = find(ismember(model.rxns, rxnsub));
        model.lb(f) = 0;
        
        % Modify the model by setting the lower bound of 'rxnslb' reactions to 0
        f = find(ismember(model.rxns, rxnslb));
        model.ub(f) = 0;

        % Run the task 'Tasks_carbon.xlsx' again to get the new essential reactions and fluxes
        [taskReport, essentialReactons, taskStructure, essentialFluxes] = checkTasks(model, 'Tasks_Co2.xlsx', true, 'getEssential', true);
    end

    % Assign the final modified model to 'modifiedModel'
    modifiedModel = model;
    
    % Restore the bounds of 'EXC_BOTH_C00007_e' in the model
    f = find(contains(model.rxns,'EXC'));
    modifiedModel.lb(f) = 0;
    modifiedModel.ub(f) = 0;
    f = find(contains(model.rxns,EXC));
    modifiedModel.lb(f) = -20;
    modifiedModel.ub(f) = 50;

    % Merge the lists of reaction IDs, formulas, and fluxes
    rxnall = [rxnsub; rxnslb];
    formulasAll = [formulasSub; formulasSlb];
    fluxAll = [fluxesub; fluxeslb];
    
    % Display the total number of modified reactions
    fprintf('Number of carbon containing reactions modified: %d\n', length(rxnall))
    
    % Write the modified reactions, their formulas, and their fluxes to an Excel file
    rxnFormulasFluxes = [rxnall, formulasAll, num2cell(fluxAll)];
    T = cell2table(rxnFormulasFluxes, 'VariableNames', {'Reaction', 'Formula', 'Flux'});
    writetable(T, 'carbon containing reactions modified.xlsx');
    
    % Display a message indicating that the Excel file has been created
    display("The reactions modified with their details are saved in the excel file 'carbon containing reactions modified.xlsx'")
end
