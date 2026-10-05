    modelM = model_red;
    obj = find(modelM.c);
    bbbs_m = modelM.mets(find(modelM.S(:,find(modelM.c)) < 0));
    r = length(modelM.rxns);
    modelM = addDemandReaction(modelM, bbbs_m);
    
    sol_bbb = zeros(length(bbbs_m), 1);  % Initialize sol_bbb
    
    for i = 1:length(bbbs_m)
        try
            modeli = modelM;
            modeli.c = zeros(length(modeli.rxns), 1);
            modeli.c(r + i) = 1;
            
            % Solve the model
            sol = optimizeCbModel(modeli, 'max');
            
            % Check solution value and update
            if sol.f < 10^-6
                sol.f = 0;
            end
            
            % Store the solution flux value for this metabolite
            sol_bbb(i, 1) = sol.f;
            
        catch ME
            % If an error occurs, print a warning and continue with the next iteration
            warning('Error with metabolite %d: %s. Skipping this iteration.', i, ME.message);
            sol_bbb(i, 1) = 0;  % Assign zero if the optimization fails
        end
    end