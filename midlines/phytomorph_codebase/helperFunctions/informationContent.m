function [c,C] = informationContent(target,source)
    
    TN = single((target == 0))'*single((source == 0));
    FP = single((target == 0))'*single((source == 1));
    FN = single((target == 1))'*single((source == 0));
    TP = single((target == 1))'*single((source == 1));
   

    C(1,1) = TN;
    C(1,2) = FN;
    C(2,1) = FP;
    C(2,2) = TP;
    
    TOP = prod(diag(C)) - prod(diag(imrotate(C,90)));
    BOTTOM = prod([sum(C,1) sum(C,2)']).^.5;
    c = TOP.*BOTTOM.^-1;
    
end