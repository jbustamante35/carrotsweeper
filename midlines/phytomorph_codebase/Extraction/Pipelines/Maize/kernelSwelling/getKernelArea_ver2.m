function [boundary currentArea F S] = getKernelArea_ver2(I,areaThresh,dB,F,S,fr,thresh)   
    dispi = 0;
    try
        if fr <= thresh
            I = imfilter(I,fspecial('disk',7),'replicate');
            R = rgb2hsv_fast(I,'single','H');
            R = rem(R + .4,1);
            level = graythresh(R);        
            MASK = R < mean(level);
            E = strel('disk',7);
            MASK = bwareaopen(MASK,areaThresh);
            MASK = imclearborder(MASK);
            MASK = imclose(MASK,E);
            MASK = imfill(MASK,'holes');
            B = bwboundaries(MASK);
           
            
          
        else
            lambda = linspace(-5,5,10);
            PDF_inner = gmdistribution.fit(squeeze(S(:,:,1)),1);
            PDF_outer = gmdistribution.fit(squeeze(S(:,:,2)),1);
            for e = 1:size(dB,1)
                for l = 1:numel(lambda)
                    ptINNER(e,l,:) = dB(e,:) + F(e,:,2)*(lambda(l)-5);
                    ptOUTER(e,l,:) = dB(e,:) + F(e,:,2)*(lambda(l)+5);
                end
            end
            for k = 1:size(I,3)
                valueINNER(:,:,k) = ba_interp2(double(I(:,:,k)),ptINNER(:,:,2),ptINNER(:,:,1));
                valueOUTER(:,:,k) = ba_interp2(double(I(:,:,k)),ptOUTER(:,:,2),ptOUTER(:,:,1));
            end
            SZ = size(valueINNER);
            probINNER = pdf(PDF_inner,reshape(valueINNER,[SZ(1)*SZ(2) SZ(3)]));
            probOUTER = pdf(PDF_outer,reshape(valueOUTER,[SZ(1)*SZ(2) SZ(3)]));
            probINNER = reshape(probINNER,SZ(1:2));
            probOUTER = reshape(probOUTER,SZ(1:2));
            PROB = log(probINNER) + log(probOUTER);
            [V,sidx] = max(PROB,[],2);
            sidx = round(imfilter(sidx,ones(5,1)/5,'replicate'));
            % build new contour
            for e = 1:size(dB,1)
                ndB(e,:) = dB(e,:) + lambda(sidx(e))*F(e,:,2);
            end
            %ndB = round(ndB);
            %ndB = unique(ndB,'rows');
            MASK = poly2mask(ndB(:,2),ndB(:,1),size(I,1),size(I,2));
            %MASK = imfilter(MASK,fspecial('average',[5 5]));
            MASK = bwareaopen(MASK,100);
            B = bwboundaries(MASK);
            %B = {ndB};
        end
        currentArea = sum(MASK(:));
        boundary = B{1};
        % create frame, sample frame
        F = calcFrame(B{1},3);
        if fr <= thresh
            tmpS = sampleFrame(B{1},F,double(I),5);
            S = cat(1,S,tmpS);
        end
        if size(boundary,1) < 10
            here= 1;
        end
        if dispi
            imshow(MASK,[]);
            drawnow
        end
    catch ME
        boundary = [];
        currentArea = 0;
    end
end

function [F] = calcFrame(dB,para)
        h0 = fspecial('gaussian',[5*para 1],para);
        h1 = gradient(h0);
        h2 = gradient(h1);
        h1 = h1 / sum(h1);
        h2 = h2 / sum(h2);
        % calculate curvature            
        tmpT = imfilter(dB,h1,'circular');
        L = sum(tmpT.*tmpT,2).^.5;
        tmpT = bsxfun(@times,tmpT,L.^-1);
        tmpN = [-tmpT(:,2),tmpT(:,1)];
        F = cat(3,tmpT,tmpN);
end

function [S] = sampleFrame(dB,F,I,lambda)
    for e = 1:size(dB,1)
        ptL(e,:) = dB(e,:) - lambda*F(e,:,2);
        ptG(e,:) = dB(e,:) + lambda*F(e,:,2);
    end
    for k = 1:size(I,3)
        S(:,k,1) = ba_interp2(I(:,:,k),ptL(:,2),ptL(:,1));
        S(:,k,2) = ba_interp2(I(:,:,k),ptG(:,2),ptG(:,1));
    end
end