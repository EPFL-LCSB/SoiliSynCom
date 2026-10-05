% Define the function to modify the model
function model = EXCTransAdder(model, EXC, EXCFormula, Trans, TransFormula)
    stoich = extract_stoichiometries(TransFormula);
    for i = 1:length(EXC)
        model = addReaction(model, EXC{i},'reactionFormula',EXCFormula{i,1},...
            'reversible',1,'lowerBound',0,'upperBound',0,...
            'subSystem','ExchangeReaction');
    end
    for i = 1:length(Trans)
        model = addReaction(model,Trans{i},'reactionFormula',TransFormula{i,1},...
            'reversible',1,'lowerBound',-50,'upperBound',50,...
            'subSystem','Transport','stoichCoeffList',stoich{i});
    end
    f = find(contains(model.mets, '_c'));
    model.metComps(f,:) = 1;
    model.metCompSymbol(f,:) = {'c'};
    f = find(contains(model.mets, '_p'));
    model.metComps(f,:) = 2;
    model.metCompSymbol(f,:) = {'p'};
    f = find(contains(model.mets, '_e'));
    model.metComps(f,:) = 3;
    model.metCompSymbol(f,:) = {'e'};
    model.metSEEDID = model.mets; 
    model.metSEEDID = strrep(model.metSEEDID,'_c',''); 
    model.metSEEDID = strrep(model.metSEEDID,'_e',''); 
    model.metSEEDID = strrep(model.metSEEDID,'_p',''); 
    model.metSEEDID = strrep(model.metSEEDID,' ','');
end
