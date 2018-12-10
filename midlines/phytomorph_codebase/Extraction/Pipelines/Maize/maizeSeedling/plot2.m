function [] = plot2(M1,M2,G,Genotype,T,Control,color)
    hold on
    g1 = strcmp(G,Genotype) & strcmp(T,Control);
    d1 = mean(M1(:,g1),2);
    d2 = mean(M2(:,g1),2);
    s1 = std(M1(:,g1),1,2)*size(M1,2)^-.5;
    s2 = std(M2(:,g1),1,2)*size(M2,2)^-.5;
    errorbar(d1,s1,color);
end