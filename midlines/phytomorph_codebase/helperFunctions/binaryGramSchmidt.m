function [G] = binaryGramSchmidt(basis,vector)
    % positive gain
    G(:,1) = basis==0 & vector==1;
    % negative gain
    G(:,2) = basis==1 & vector==0;
end