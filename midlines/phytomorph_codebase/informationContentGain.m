function [cG,C,totC] = informationContentGain(target,source,basis)

    [c,originalC] = informationContent(target,source);

    %[G] = binaryGramSchmidt(basis,source);
    G = all(bsxfun(@and,basis,source),2);
    
    
    
    
    TN = single((target == 0))'*single((G == 0));
    FP = single((target == 0))'*single((G == 1));
    FN = single((target == 1))'*single((G == 0));
    TP = single((target == 1))'*single((G == 1));
   

    C(1,1) = TN;
    C(1,2) = FN;
    C(2,1) = FP;
    C(2,2) = TP;
    %C = originalC + C;
    
    
    TOP = prod(diag(C)) - prod(diag(imrotate(C,90)));
    BOTTOM = prod([sum(C,1) sum(C,2)']).^.5;
    cG = TOP.*BOTTOM.^-1;
    totC = cG;


    cG = cG - c;
    
end