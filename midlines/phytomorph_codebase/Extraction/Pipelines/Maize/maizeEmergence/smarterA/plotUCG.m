function [] = plotUCG(pt,edge)
    plot(pt(:,1),pt(:,2),'k.');
    for e = 1:size(edge,3)
        plot(edge(:,1,e),edge(:,2,e),'k');
    end
end