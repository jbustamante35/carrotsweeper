function [grade,ret1,ret2,BO,TC] = opti(source,target,para,oI,init,delta,confusionFunction)
    if nargin <= 4 | (isempty(init) & isempty(delta))
        init = [.1 1*3/9 1*3/9 7*3/9 100 .3 30 3 3];
        delta = [.05 .2 .2 .2 10 0.15 10 1 1]*.05^-1;
    end

    TIN = [];
    TC = 0;
    grade = NaN;
    ret1 = NaN;
    %para = para.*alpha.^-1;
   
    para = init + (para - delta);
    ag = [];
    ret1 = [];
    BO = [];
    disp = false;
    %para(5) = 0;
    %targetSZ = size(target);
    %TOTX = zeros(prod(targetSZ(1:2)),1);
    %TOTY = zeros(prod(targetSZ(1:2)),1);
    %str = 1;
    %skip = prod(targetSZ(1:2));
    %stp = str + skip - 1;
    
    parfor o = 1:size(source,4)
        %str = (o-1)*skip + 1;
        %stp = str + skip - 1;
        tmpS = source(:,:,:,o);
        %tmpS(:,:,4) = std(tmpS,1,3);
        
        for slice = 1:size(tmpS,3)
            tmpS(:,:,slice) = tmpS(:,:,slice) * para(slice+1);
        end
        
        tmpS = imfilter(tmpS,fspecial('gaussian',[21 21],round(para(8))),'replicate');
        
        
        % try giving it a chance by two
        %tmpS = sort(tmpS,3,'descend');
        
        
        marker = mean(tmpS,3);
       
        
        
        LEVEL = para(1);
        mm = marker - LEVEL;
        RECON = imreconstruct(mm,marker);
        BO = (mm - RECON) > - LEVEL;
        
        BO = imfill(BO,'holes');
        
        %{
        BBO = imclearborder(BO) == 0 & BO == 1;
        BO = bwareaopen(BO,round(para(5)));
        BBO = bwareaopen(BBO,round(para(7)));
        BO = logical(BO + BBO);
        %}
        
        BO = imclose(BO,strel('disk',5,0));
        %BO = bwareaopen(BO,round(para(5)));
        
        POT = double(BO).*marker;
        
        BO = POT > para(6);
        POT = POT.*BO;
        di = round(para(9));
        if di < 0
            di = 0;
        end
        BO = imdilate(BO,strel('disk',di,0));
        
        
        BBO = imclearborder(BO) == 0 & BO == 1;
        BO = bwareaopen(BO,round(para(5)));
        BBO = bwareaopen(BBO,round(para(7)));
        BO = logical(BO + BBO);
        
        R = regionprops(BO);
        N(o) = numel(R);
        
        
        
        if ~isempty(target)
            tmpT = logical(target(:,:,o));
            tmpT = imerode(tmpT,strel('disk',2,0));
        end
        
        RR = regionprops(logical(BO),'Area','MajorAxisLength');
        
        %c1 = count([RR.Area]);
        %c2 = count([RR.MajorAxisLength]);
        %IDX2 = c1 == 2 | c2 == 2;
        %TC = sum(c1);
        
        out = flattenMaskOverlay(oI(:,:,o),BO,.5,'r');
        
        
        if ~isempty(target)
            out = flattenMaskOverlay(out,tmpT,.5,'b');
        end
        
        
        ret2(:,:,:,o) = out;
        
        
        if o == 1
            if  disp
                imshow(out,[]);
                title(num2str(o));
                drawnow
            end
           
        end
        
        if ~isempty(target)
            
            
            
            
            tmp = target(:,:,o);
            tmp = tmpT;
            
            %TIN = [TIN;tmpT(:)];
            %TOUT = [TOUT;POT(:)];
            
            targets = [tmp(:) ~tmp(:)];
            outputs = [BO(:) ~BO(:)];
            
            HOPE = marker.*BO;
            TOTX{o} = targets(:,1);
            TOTY{o} = HOPE(:);
            %[tpr,fpr] = roc(double(targets(:,1))',double(outputs)');
            %ag = [ag;[tpr(2) fpr(2)]];
            
            
            %C = confusionmat(logical(targets(:,1)),logical(outputs(:,1)),'order',logical([0 1]));
            gr2(o) = confusionFunction(targets,outputs);
            %gr2(o) = TOP*BOTTOM^-1;
            
            %gr(o) = targets(:,1)'*outputs(:,1);
%            grade(o) = tpr(2).*fpr(2).^-1;


            %ret1 = [ret1;[tpr(2),fpr(2)]];
        end
        
        
    end
    
    if ~isempty(target)
        
        %[Xp,Yp] = perfcurve(TOTX,TOTY,1);
        if N(1) == 1
            here = 1;
        end
        imshow(ret2(:,:,:,1),[]);
        drawnow
        %waitforbuttonpress
        %-(mean(grade))
        %grade = -(mean(grade));
        % grade;
        rmidx = isinf(gr2);
        gr2(rmidx) = NaN;
        MX = max(gr2);
        gr2(isnan(gr2)) = MX;
        %grade = -mean(gr);
        grade = -mean(gr2);
        if grade < -100
            here = 1;
        end
    end
   
    
    
    
    
    
    
        
        
end