function [] = plotX(M1,G,Genotype,T,Control,color)
    hold on
    g1 = strcmp(G,Genotype) & strcmp(T,Control);
    d1 = mean(M1(:,g1),2);
    s1 = std(M1(:,g1),1,2)*size(M1,2)^-.5;
    errorbar(d1,s1,color);
end