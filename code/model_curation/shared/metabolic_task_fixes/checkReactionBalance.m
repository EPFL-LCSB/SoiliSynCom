function isBalanced = checkReactionBalance(S, elementMatrix, reactionIndex)
    % S - stoichiometric matrix (m x n, where m is the number of metabolites and n is the number of reactions)
    % elementMatrix - element matrix (m x p, where m is the number of metabolites and p is the number of elements)
    % reactionIndex - index of the reaction to check

    % Get the stoichiometric vector for the reaction
    stoichiometricVector = S(:, reactionIndex);

    % Compute the element balance for the reaction
    elementBalance = elementMatrix' * stoichiometricVector;

    % Check if the reaction is balanced (i.e., the element balance vector is all zeros)
    isBalanced = all(elementBalance == 0);
end
