function [patchStructure,pTable] = extractPixelListsFromPatches(patchStructure,imageSize,fileName)
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['Start pixel extration\n']);tm = clock;
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['start:Start patch iteration.\n']);
    
    tm1 = [];
    tm2 = [];
    pTable = table();
    cnt = 1;
    
    % extract features from each potential spore
    for e = 1:numel(patchStructure)
        try
            % get the patch size
            patchSZ = (size(patchStructure(e).patch)-1)/2;
            % get the patch location
            tmpLoc = flipdim(patchStructure(e).patchLocation,2) - patchSZ;
            % get the patch center
            patchCenter = (size(patchStructure(e).mask)-1)/2;
            
            
            
            %fprintf(['start:Calculating spore pixels.\n']);tic;
            tic;
            % fill holes in patch mask
            tmpM = imfill(patchStructure(e).mask,'holes');
            % find holes
            msk = imfill(~tmpM,patchCenter+1) == 1 & patchStructure(e).mask == 1;
            midx = find(msk);
            
            
            subIndex = [];
            [subIndex(:,1),subIndex(:,2)] = ind2sub(size(tmpM),midx);
            subIndex = bsxfun(@plus,subIndex,tmpLoc);
            subIndex = sub2ind(imageSize,subIndex(:,1),subIndex(:,2));
            patchStructure(e).sporeIDX = subIndex;
            patchStructure(e).pureSporeMask = msk;
            
            
            % compute the features of the patch mask
            sporeProps = regionprops(logical(patchStructure(e).pureSporeMask),'PixelIdxList','Eccentricity','MajorAxisLength','MinorAxisLength');
            %fprintf(['end:Calculating spore pixels.:' num2str(toc) '\n']);
            if (mod(e,10)==0);fprintf(['.\n']);else;fprintf(['.']);end
            tm1 = [tm1 toc];


            %fprintf(['start:Calculating surface kurvature.\n']);tic;
            tic;
            para.scales.value = 2;
            para.resize.value = 1;
            [K] = surKur(patchStructure(e).patch,para);
            K = bindVec(K(:,:,1));
            %fprintf(['end:Calculating surface kurvature.:' num2str(toc) '\n']);
            if (mod(e,10)==0);fprintf(['.']);else;fprintf(['.']);end
            tm2 = [tm2 toc];

            
            

            K(midx) = mean(K(:));
            T = adaptthresh(K,.4);

            hyphaMask = imbinarize(K,T);
            Rh = regionprops(hyphaMask,'PixelIdxList');
            Rs = regionprops(imdilate(msk,strel('disk',3,0)),'PixelIdxList');
            H = zeros(size(hyphaMask));
            hyphaIDX = [];
            for h = 1:numel(Rh)
                if ~isempty(intersect(Rs(1).PixelIdxList,Rh(h).PixelIdxList))
                    H(Rh(h).PixelIdxList) = 1;

                    subIndex = [];
                    [subIndex(:,1) subIndex(:,2)] = ind2sub(size(tmpM),Rh(h).PixelIdxList);
                    subIndex = bsxfun(@plus,subIndex,tmpLoc);
                    subIndex = sub2ind(imageSize,subIndex(:,1),subIndex(:,2));
                    hyphaIDX = [hyphaIDX;subIndex];
                end
            end
            patchStructure(e).hyphaMask = hyphaMask;
            patchStructure(e).purehyphaMask = logical(H);
            patchStructure(e).hyphaIDX = subIndex;


            totalHypaArea = sum(patchStructure(e).purehyphaMask(:));
            totalSporeArea = sum(patchStructure(e).pureSporeMask(:));
            DT = bwdist(patchStructure(e).pureSporeMask);
            hidx = find(patchStructure(e).purehyphaMask);
            meanDistance = mean(DT(hidx));
            stdDistance = std(DT(hidx),1,1);


            
            
            pTable{cnt,'fileName'} = {fileName};
            pTable{cnt,'totalHypaArea'} = totalHypaArea;
            pTable{cnt,'totalSporeArea'} = totalSporeArea;
            pTable{cnt,'meanHyphaDistance'} = meanDistance;
            pTable{cnt,'stdHyphaDistance'} = stdDistance;
            pTable{cnt,'rowIndex'} = patchStructure(e).patchLocation(2);
            pTable{cnt,'columnIndex'} = patchStructure(e).patchLocation(1);
            pTable{cnt,'species'} = {patchStructure(e).species};
            pTable{cnt,'sporeEccentricity'} = sporeProps(1).Eccentricity;
            pTable{cnt,'sporeMajorAxisLength'} = sporeProps(1).MajorAxisLength;
            pTable{cnt,'sporeMinorAxisLength'} = sporeProps(1).MinorAxisLength;
            % make the feature vector based on distance from outline
            dVec = linspace(0,30,31);
            for d = 1:(numel(dVec)-1)
                dMSK = DT > 0 & DT >= dVec(d) & DT <= dVec(d+1);
                fidx = find(dMSK);
                pTable{cnt,['K' num2str(dVec(d)) '_' num2str(dVec(d+1))]} = mean(K(fidx));
            end

            
            cnt = cnt + 1;

        catch ME
            what = 1;
        end
    end
    fprintf(['\nreport:Calculating spore pixels time:' num2str(mean(tm1)) '\n']);
    fprintf(['report:Calculating surface kurvature:' num2str(mean(tm2)) '\n']);
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
    fprintf(['End pixel extration.' num2str(etime(clock,tm)) '\n']);
    fprintf(['*******************************************\n']);
    fprintf(['*******************************************\n']);
end