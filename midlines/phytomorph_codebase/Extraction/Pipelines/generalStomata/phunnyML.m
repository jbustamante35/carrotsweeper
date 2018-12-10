function [d] = phunnyML(I,P,D,p,z,U,E)


    [T] = generateDomainTransformation(p(1:5));
    [D] = transformDomain(D,T);
    [sI] = myInterp2Sampler(I,P,D,z);
    
    C = PCA_REPROJ_T(sI(:),E,U);
    ssI = PCA_BKPROJ_T(C,E,U);
    d = norm(ssI(:) - sI(:));
    
    
end