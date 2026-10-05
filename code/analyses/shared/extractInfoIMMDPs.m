function [NumActive, StatsMets, Mets] = extractInfoIMMDPs(model, DPs)
% Extract statistics of metabolites in DPs of IMM/IMS
%
% USAGE:
%
%    [NumActive, StatsMets, Mets] = extractInfoIMMDPs(model, DPs)
%
% INPUT:
%    model:           TFA model with MILP structure for IMM analysis
%    DPs:             Directionality profile matrix with alternatives in
%                     each column from IMM analysis 
%
% OUTPUTS:
%    NumActive:       Number of mets in the IMM of each alternative
%    StatsMets:       Statistics on appearence of mets among all
%                     alternatives
%    Mets:            Summary of drains and metabolites
%
% .. Author:
% Anush Chiappino-Pepe 2014
% 

drains = getAllVar(model,{'BFUSE'});
NumActive = zeros(1,size(DPs,2));
StatsMets = zeros(length(drains),1);

% extract from DPs the matrix with the information for the drains
MatrixInfo = DPs(size(DPs,1)-length(drains)+1:size(DPs,1),:);

% identify number of mets in the IMM of each alternative
for i = 1:size(MatrixInfo,2)
    NumActive(i) = length(drains)-sum(MatrixInfo(:,i));
end

% get statistics on appearence of mets among all alternatives
for i = 1:size(MatrixInfo,1)
    StatsMets(i) = size(DPs,2)-sum(MatrixInfo(i,:));
end

% get met names associated to statistics (should be also equivalent to the
% output drains of analysisIMM.m)
rxndrains = model.varNames(drains);
rxndrains = strrep(rxndrains,'BFUSE_','');
rxndrains = strrep(rxndrains,'R_','');
rxndrains = strrep(rxndrains,'F_','');
Mets(:,1) = printRxnFormula(model,rxndrains);
Mets(:,2) = printRxnFormula(model,rxndrains,0,0,1);
end