function [X] = extractionAtom0(X,alongDim,numberFreq)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alongDim :=  operate along dim
    % numberFreq := number of frequences to extract
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    M = mfftm(size(X,alongDim),numberFreq);
    if alongDim == 1
        X = mtimesx(M,X,'T');
    else
        X = mtimesx(M,X);
    end
end