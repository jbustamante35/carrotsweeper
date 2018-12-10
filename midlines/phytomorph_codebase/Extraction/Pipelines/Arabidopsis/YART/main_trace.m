function [pheno] = main_trace(I,X,disp,toC)
    
    pheno = myTraceLow(I,X,disp,toC);
    para.SNIP = 10;
    pheno = measurePheno(pheno,para);
    phenoView(I,pheno);
end