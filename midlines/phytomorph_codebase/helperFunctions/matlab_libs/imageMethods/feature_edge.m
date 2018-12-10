function [X] = feature_edge(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% input para
    E = edge(I,para{1});
    [r c] = find(E);
    X = [r c];
end