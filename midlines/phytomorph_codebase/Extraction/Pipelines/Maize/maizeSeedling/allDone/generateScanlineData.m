function [tmp_HLINE_Y,tmp_HLINE_X,tmp_VLINE_Y,tmp_VLINE_X] = generateScanlineData(I,tmpMask,scanlineSZ,imgSZ)
        tmp_HLINE_Y = [];
        tmp_HLINE_X = [];
        tmp_VLINE_Y = [];
        tmp_VLINE_X = [];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(tmpMask)
            TMP = [];
            for k = 1:size(tmpMask,3)
                tmp = im2colF(tmpMask,[imgSZ(1) scanlineSZ],[1 1]);
                tmp = reshape(tmp,[imgSZ(1) scanlineSZ 1 size(tmp,2)]);
                tmp_HLINE_Y = [tmp_HLINE_Y any(squeeze(tmp(:,(end-1)/2,1,:)),1)'];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(I)
            TMP = [];
            for k = 1:3
                tmp = im2colF(I(:,:,k),[imgSZ(1) scanlineSZ],[1 1]);
                tmp = reshape(tmp,[imgSZ(1) scanlineSZ 1 size(tmp,2)]);
                TMP = cat(3,TMP,tmp);
            end
            tmp_HLINE_X = TMP;
            I = permute(I,[2 1 3]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(tmpMask)
            tmpMask = permute(tmpMask,[2 1]);
            TMP = [];
            for k = 1:size(tmpMask,3)
                tmp = im2colF(tmpMask,[imgSZ(2) scanlineSZ],[1 1]);
                tmp = reshape(tmp,[imgSZ(2) scanlineSZ 1 size(tmp,2)]);
                tmp_VLINE_Y =[tmp_VLINE_Y any(squeeze(tmp(:,(end-1)/2,1,:)),1)'];
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(I)
            TMP = [];
            for k = 1:3
                tmp = im2colF(I(:,:,k),[imgSZ(2) scanlineSZ],[1 1]);
                tmp = reshape(tmp,[imgSZ(2) scanlineSZ 1 size(tmp,2)]);
                TMP = cat(3,TMP,tmp);
            end
            tmp_VLINE_X = TMP;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end