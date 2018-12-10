function [fStack,subIR] = myCropper(I,BOX,scale,verbose,useScale)
    if nargin == 3;verbose = false;end
    
    try
        fStack = {};
        for e = 1:size(BOX,1)
            if verbose;tic;end
            
            
            
            
            subI = mimcrop(I,BOX,round(scale*BOX));
            
            
            
            %subI = imcrop(I,BOX(e,:));
            %subI = I(round(BOX(e,2):(BOX(e,2)+BOX(e,4))),round(BOX(e,1):(BOX(e,1)+BOX(e,3))),:);
            %subI = I(round(BOX(e,2):((BOX(e,2)+BOX(e,4)-1))),round(BOX(e,1):((BOX(e,1)+BOX(e,3)-1))),:);
            
            %{
            if nargin == 5  
                subI = imresize(subI,fliplr(useScale));
            end
            %}
            
            
            tStack = featureFunction_wt(subI,e,scale,3,'haar',false);
            
            
            
            
            
            
            
            %tStack = featureFunction(subI,e);

            if e == 1
                for n = 1:numel(tStack)
                    fStack{n} = [];
                end
            end

            for n = 1:numel(tStack)
                fStack{n} = cat(4,fStack{n},tStack{n});
            end
            if verbose;fprintf(['done with:' num2str(e) ':' num2str(size(BOX,1)) '@' num2str(toc) '\n']);end
            subIR(:,:,:,e) = subI;
        end
    catch ME
        ME
    end
end