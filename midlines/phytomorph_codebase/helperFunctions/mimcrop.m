function [subI,rotM,U] = mimcrop(I,BOX,newSize,rotPot)
    rotM = 1;
    U = 0;
    if newSize ~= Inf
        
        % if new size is a scaling factor
        if numel(newSize) == 1
            newSize = BOX(3:4)*newSize;
        end
        [X(:,:,2),X(:,:,1)] = ndgrid(linspace(BOX(2),BOX(2)+BOX(4),newSize(2)),linspace(BOX(1),BOX(1)+BOX(3),newSize(1)));
        
        % place rot here
        if nargin == 4
            if ~isempty(rotPot)
                if numel(rotPot) ~= 1
                    TH = linspace(rotPot(1),rotPot(2),rotPot(3));
                    IDX = randi(rotPot(3),1);
                    TH = TH(IDX);
                else
                    TH = rotPot;
                end
                
                rotM = [[cos(TH) sin(TH)];[-sin(TH) cos(TH)]];
                sz = size(X);
                X = reshape(X,[prod(sz(1:2)) sz(3)]);
                U = mean(X,1);
                X = bsxfun(@minus,X,U);
                X = (rotM*X')';
                X = bsxfun(@plus,X,U);
                X = reshape(X,sz);
            
            else
                sz = size(X);
                X = reshape(X,[prod(sz(1:2)) sz(3)]);
                U = mean(X,1);
                X = reshape(X,sz);
            end
        else
            sz = size(X);
            X = reshape(X,[prod(sz(1:2)) sz(3)]);
            U = mean(X,1);
            X = reshape(X,sz);
        end
        
        
        subI = ba_interp2(I,X(:,:,1),X(:,:,2));
        
    else
    
        subI = I(BOX(1):(BOX(1)+BOX(3)),BOX(2):(BOX(2)+BOX(4)),:);
        
        
    end
    
end