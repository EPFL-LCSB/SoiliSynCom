function [InorganicMetExchange,InorganicMetExchangeIdx] = findInorganicsIMEs(model,NoFormulaStr,ExtraFormulas,ToExclude)
% Simple function to find inorganic metabolites from COBRA compatible
% models that have a metFormula entry (model.metFormulas).
%
% INPUTS
% - model: A COBRA compatible model structire
% - NoFormulaStr: A cell array containing characteristic strings used in
%   model.metFormulas to indicate that no chemical formula is available for
%   a particular metabolite. This is model specific (default is: 'NA',
%   'None','NULL','').
% - ExtraFormulas: To the list of inorganic compounds, we might also want
%   to add some compounds that are not captured by the simplistic approach
%   of this function. If not set, by default we also add cyanate (CNO),
%   bicarbonate (CHO3), carbon dioxide (CO2), thiocyanate (CNS), and
%   hydrogen cyanide (CHN). 
% - FormulasToExclude: some metabolites whose formula is misleading, and might end
%   up in the inorganic list, such as:
%   S2Fe2: [2Fe-2S] iron-sulfur cluster, which is not inorganic
% 
% OUTPUTS
% - InorgMetSEEDIDs: SEEDIDs of all the inorganic metabolites with no
%   repetitions.
% - InorgMets: cell array of size n-by-5, where n is the number of
%   inorganic compounds found in the model. The 5 columns correspond to:
%   1) Index-number of metabolites in the model metabolite list
%   2).mets
%   3).metNames
%   4).metFormulas
%   5).metSEEDID
% 
% LCSB, EPFL, December 2014


if nargin<4
    FormulasToExclude = {'S2Fe2','S4Fe4'};
end

if nargin<3
    ExtraFormulas = {'CNO','CHO3','CO2','CNS','CHN','H2S2X','S2X','CHO'};
end

if nargin<2
    NoFormulaStr = {'NA','None','NULL',''};
end
    
% Make sure that the provided model contains a field with the chemical formulas
model_fields = fieldnames(model);

if ~ismember(model_fields,'metFormulas')
    error('Model File Does not contain Metabolite Formulas or the variable is not named model.metFormulas')
end

% Find the indices of the metabolites whose chemical formulas does not
% contain 'C' followed by another capital letter or any number. For example
% indices of elements such as 'Ca', 'Cd', etc will be found:
InorganicMetID = find(cellfun(@isempty,regexp(model.metFormulas,'(?-i)C[A-Z_0-9]')));
% Attention: are there any falsely identified organic metabolites?
FalseInorgMetIds = intersect(InorganicMetID,find_cell(FormulasToExclude,model.metFormulas));
if ~isempty(FalseInorgMetIds)
    InorganicMetID = setdiff(InorganicMetID,FalseInorgMetIds);
end
% From these indices (InorganicMetID) exclude the indices that correspond
% to "NoFormulaStr".
NoFormulaStrID = find(ismember(model.metFormulas,NoFormulaStr));
InorganicMetID = setdiff(InorganicMetID,NoFormulaStrID);
% To these indices (InorganicMetID) add the indices corresponding to the 
% "ExtraFormulas"
ExtraFormulasID = find(ismember(model.metFormulas,ExtraFormulas));
InorganicMetID = union(InorganicMetID,ExtraFormulasID);

if size(InorganicMetID,2) > size(InorganicMetID,1)
    InorganicMetID = InorganicMetID';
end

mets = model.mets(InorganicMetID);
mets = mets(find(contains(mets,'_e')));

[~,InorganicMetExchangeIdx] = ismember(strcat('EXC_BOTH_',mets),model.rxns);
InorganicMetExchangeIdx = InorganicMetExchangeIdx(find(InorganicMetExchangeIdx));
InorganicMetExchange = model.rxns(InorganicMetExchangeIdx);

end