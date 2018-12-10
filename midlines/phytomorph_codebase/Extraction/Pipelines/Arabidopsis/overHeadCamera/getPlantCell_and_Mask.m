function [Area,QR_text,subI,plantMask_tmp] = getPlantCell_and_Mask(I,box,GMModel,map)
    for b = 1:numel(box)

        %subI_tmp = imcrop(I,box{b});
        subI_tmp = mimcrop(I,box{b},1,[]);
        sz = size(subI_tmp);
        cidx = GMModel.cluster(reshape(subI_tmp,[prod(sz(1:2)) sz(3)]));
        subL = reshape(cidx,sz(1:2));

        try
            label_nodither = rgb2ind(subI_tmp,map,'nodither');
            white = label_nodither == 3;
            white = bwareaopen(white,1000);
            white = bwlarge(white);
            Rwhite = regionprops(white);
            label = imcrop(subI_tmp,Rwhite(1).BoundingBox);
            QR_text{b} = decode_qr(label);
        catch ME
            getReport(ME)
            QR_text{b} = '';
        end



        plantMask_tmp = subL == 3;
        plantMask_tmp = bwareaopen(plantMask_tmp,300);
        plantMask_tmp = imclearborder(plantMask_tmp);

        Area(b) = sum(plantMask_tmp(:));

        if nargout > 1
            subI(:,:,:,b) = subI_tmp;
            plantMask(:,:,b) = plantMask_tmp;
        end
    end
end