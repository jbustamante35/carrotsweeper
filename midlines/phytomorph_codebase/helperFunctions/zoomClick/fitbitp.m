function [Y] = fitbitp(X,bx,map,beta,WINDOW,raster)
    if nargin == 5
        raster = false;
    end
    for tr = 1:size(X,4)
        %dX(:,:,tr) = dither(X(:,:,:,tr), map);
        dX(:,:,tr) = rgb2ind(X(:,:,:,tr), map,'nodither');
    end
    
    [bX] = bits(dX,bx);
    
    if raster
        dX = double(dX);
        B = im2colF((dX(:,:,1)),[15 15],[1 1]);
        B = zeros(size(B,1),size(B,2),size(dX,3));
        for tr = 1:size(dX,3)
            %tic;
            B(:,:,tr) = im2colF((dX(:,:,tr)),[15 15],[1 1]);
            %tm(tr) = toc;
            %mean(tr)*(size(dX,3)-tr)/60
            tr
        end
        B = squeeze(sum(B,1));
        bX = B;
    end
    
    
    
    
    abX = [ones(size(bX,1),1) bX];
    %{
    OV = zeros([size(abX,1) size(beta,2) size(beta,3)]);
    for e = 1:size(beta,3)
        OV(:,:,e) = mtimesx(abX,beta(:,:,e));  
    end
    %}
    OV = mtimesx(abX,beta);
    
    for tr = 1:size(OV,1)
        for d = 1:size(OV,3)
            fidx = find(OV(tr,:,d) > .5);
            if isempty(fidx)
                fidx = 1;
            end
            Y(tr,d) = WINDOW(fidx(end),d);
        end
    end

end