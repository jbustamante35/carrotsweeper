function [B] = clearAllBorder(B)    
    G = imclearborder(B);
    B = G ~= B;
end