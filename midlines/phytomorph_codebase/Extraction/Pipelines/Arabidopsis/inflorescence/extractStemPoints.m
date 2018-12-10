function [idxL] = extractStemPoints(fileList,disp)

    para.scales.value = 2;
    para.resize.value = .75;
    idxL = [];
    for e = 1:numel(fileList)
        
        
        I = double(imread(fileList{e}))/255;
        
        
        M(:,:,e) = I;
        fM(:,:,e) = imfilter(I,fspecial('disk',21),'replicate');
        BK = imfilter(I,fspecial('disk',101),'replicate');
        sM(:,:,e) = I - BK;




        K = surKur(sM(:,:,e),para);
        sig1 = bindVec(K(:,:,1));
        stem1 = sig1 > graythresh(sig1);

        sig2 = bindVec(sM(:,:,e));
        stem2 = sig2 < graythresh(sig2);


        stemT = logical(stem1.*stem2);
        
        %stemT = imclose(stemT,strel('disk',5,0));
        
        R = regionprops(stemT,'PixelIdxList','MajorAxisLength','MinorAxisLength','Orientation');
        RA = [R.MajorAxisLength].*[R.MinorAxisLength].^-1;
        fidx = [];

        if e == 1
            fidx = find(abs([R.Orientation]) < 45 & RA > 5 & [R.MajorAxisLength] > 100);
            [~,didx] = sort(RA(fidx),'descend');
            for n = 1:numel(fidx)
                idxL = [idxL ;[R(fidx(n)).PixelIdxList e*ones(numel(R(fidx(n)).PixelIdxList),1)]];
            end
        else
            for n = 1:numel(R)
                if ~isempty(intersect(R(n).PixelIdxList,idxL(:,1)))
                    idxL = [idxL ;[R(n).PixelIdxList e*ones(numel(R(n).PixelIdxList),1)]];
                end
            end
        end

        stemT = zeros(size(stemT));
        stemT(idxL(idxL(:,2)==e,1)) = 1;

        stem(:,:,e) = stemT;
        if disp
            out = flattenMaskOverlay(I,logical(stemT));
            imshow(out,[]);
            drawnow
        end
    end
end