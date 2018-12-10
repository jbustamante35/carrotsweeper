function [fF] = featureFunction_wt(I,n,scale,L,wavelet,disp)
    TOT = 2;
    dwtmode('sp0','nodisp');
    if scale ~= 1
        I = imresize(I,scale,'nearest');
    end
    %I = ba_interp2(I,grid1(:,:,1),grid1(:,:,2));
    if disp
        imshow(I,[]);
        title(['[' num2str(size(I,1)) ',' num2str(size(I,2)) ']']);
        drawnow
    end
    [C,S] = wavedec2(I(:,:,1),L,wavelet);
    
    for l = 1:TOT
        % backward
        %[tmpH1,tmpV1,tmpD1] = detcoef2('all',C,S,TOT-(l-1));
        [tmpH1,tmpV1] = detcoef2('all',C,S,TOT-(l-1));
        % forward
        %[tmpH1,tmpV1,tmpD1] = detcoef2('all',C,S,l);
        if l ~= 1
            tmpH1 = imresize(tmpH1,[size(H,1) size(H,2)]);
            tmpV1 = imresize(tmpV1,[size(V,1) size(V,2)]);
            %tmpH1 = ba_interp2(tmpH1,grid2(:,:,1),grid2(:,:,2));
            %tmpV1 = ba_interp2(tmpV1,grid2(:,:,1),grid2(:,:,2));
            %tmpD1 = imresize(tmpD1,[size(D,1) size(D,2)]);
        end
        H(:,:,l) = tmpH1;
        V(:,:,l) = tmpV1;
        %D(:,:,l) = tmpD1;
    end
    if size(H,2) ~= 13
        here = 1;
    end
    fF{1} = H;
    fF{2} = V;
    %fF{3} = D;
end