function [tmp_squareX,tmp_squareY] = generateSquarescanData(I,tmpMask,squareSZ,PER)
    tmp_squareY = [];
    tmp_squareX = [];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(tmpMask)
        TMP = [];
        for k = 1:size(tmpMask,3)
            tmp = im2colF(tmpMask(:,:,k),[squareSZ],[1 1]);
            tmp = reshape(tmp,[squareSZ size(tmp,2)]);
            tmp = squeeze(tmp((end-1)/2,(end-1)/2,:));
            tmp_squareY = [tmp_squareY tmp(:)];
        end
        squeeze(sum(sum(tmpMask,1),2));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(I)
        TMP = [];
        for k = 1:3
            tmp = im2colF(I(:,:,k),[squareSZ],[1 1]);
            tmp = reshape(tmp,[squareSZ 1 size(tmp,2)]);
            TMP = cat(3,TMP,tmp);
        end
        tmp_squareX = TMP;
    end
    if ~isempty(PER)
        fidx1 = any(tmp_squareY==1,2);
        fidx0 = ~fidx1;
        fidx1 = find(fidx1);
        fidx0 = find(fidx0);
        ridx1 = randperm(numel(fidx1));
        ridx0 = randperm(numel(fidx0));
        n1 = round(numel(fidx1)*PER(2));
        n0 = round(numel(fidx0)*PER(1));
        tmp_squareX = cat(4,tmp_squareX(:,:,:,fidx1(ridx1(1:n1))),tmp_squareX(:,:,:,fidx0(ridx0(1:n0))));
        tmp_squareY = cat(1,tmp_squareY(fidx1(ridx1(1:n1)),:),tmp_squareY(fidx0(ridx0(1:n0)),:));
        n1
        n0
        sample = 1;
    end
    
end