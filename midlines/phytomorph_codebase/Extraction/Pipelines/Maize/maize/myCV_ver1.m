function [me] = myCV_ver1(TR,TE,n)
    trainx = TR(:,1:end-1);
    trainy = TR(:,end);
    testx = TE(:,1:end-1);
    testy = TE(:,end);
    %[BETA,SIGMA,RESID] = mvregress(trainx,trainy);
    %Ya = [testx]*BETA;
    [XL,YL,XS,YS,BETA,PCTVAR,MSE,st] = plsregress(trainx,trainy,n);
    Ya = [ones(size(testx,1),1) testx]*BETA;
    
    delta = Ya - testy;
    delta = delta.^2;
    me = [mean(delta) corr(Ya,testy)];
end