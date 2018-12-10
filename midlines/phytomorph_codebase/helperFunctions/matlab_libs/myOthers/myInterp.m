function [Y] = myInterp(I,X)
    
    Y = ba_interp2(I,X(:,1),X(:,2));
    
    
end