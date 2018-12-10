function [sidx] = findRotateFrame(G)

    %% find the rotate location
    dG = diff(G,1,3);
    dG = squeeze(sum(sum(dG.*dG,1),2));
    [~,sidx] = max(dG);
end