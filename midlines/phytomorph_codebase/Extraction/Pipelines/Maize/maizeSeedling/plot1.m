function [] = plot1(M1,M2,G,Genotype,T,Control,color)

    g1 = strcmp(G,Genotype) & strcmp(T,Control);
    d1 = mean(M1(:,g1),2);
    d2 = mean(M2(:,g1),2);
    s1 = std(M1(:,g1),1,2)*size(M1,2)^-.5;
    s2 = std(M2(:,g1),1,2)*size(M2,2)^-.5;
    plot(d1,d2,color);
    hold on
    for e = 1:numel(d1)
        E1 = linspace(-s1(e),s1(e),2) + d1(e);
        plot(E1,[d2(e) d2(e)],color);
        E1 = linspace(-s2(e),s2(e),2) + d2(e);
        plot([d1(e) d1(e)],E1,color);
    end
end