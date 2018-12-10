function [grade,myMask1,myMask2,out,ag] = generateStomataMask(maskTreeS,sMSKS,thresholdValues,oI)
    ag = [];
    for e = 1:numel(maskTreeS)
        maskTree = maskTreeS{e};
        sMSK = sMSKS(:,:,e);
        
        
        
        
        
        
        [myMask1,myMask2,totalMask] = generate12Mask(maskTree,thresholdValues,sMSK);
        
        targets = [sMSK(:)' ~sMSK(:)'];
        outputs = [totalMask(:)' ~totalMask(:)'];
        [tpr,fpr,thresholds] = roc(targets,outputs);
        ag = [ag;[tpr(2) fpr(2)]];
        grade(e) = tpr(2).*fpr(2).^-1;
        
        
        
        if nargin == 4
            out(:,:,:,e) = flattenMaskOverlay(oI(:,:,e),totalMask);
            out(:,:,:,e) = flattenMaskOverlay(out(:,:,:,e),logical(sMSK),.5,'b');
        end
        %plotroc(targets,outputs);

        %{
        TP = zeros(numel(Rtotal),1);

        for e1 = 1:numel(Rtotal)
            for e2 = 1:numel(Rclick)
                inter = intersect(Rtotal(e1).PixelIdxList,Rclick(e2).PixelIdxList);
                if ~isempty(inter)
                    TP(e1) = 1;
                end
            end
        end

        %TP:=TP==1;
        %FP:=TP==0;
        %TN:=
        %FN:=
        FN = ones(numel(Rclick),1);
        for e1 = 1:numel(Rclick)
            for e2 = 1:numel(Rtotal)
                inter = intersect(Rtotal(e2).PixelIdxList,Rclick(e1).PixelIdxList);
                if ~isempty(inter)
                    FN(e1) = 0;
                end
            end
        end








        %}


    end
    
    grade = -mean(grade);
    
    
    
    
    
    
    
    
    
end