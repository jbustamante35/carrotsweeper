function [sigX,sigY] = generateStemDomain(fileList,idxL,dataXFunc,dataYFunc,R)
    I = imread(fileList{1});
    cnt = 1;


    sigX = {};
    sigY = {};
    for e = 1:numel(fileList)
        I = double(imread(fileList{1}))/255;
        
       
        
        
        
        Z = zeros(size(I));
        fidx = find(idxL(:,2) == e);
        Z(idxL(fidx,1)) = 1;
        
        
        [D] = dataXFunc(I);
        [DM] = dataYFunc(Z);
        
        
        
        
        for r = 1:4
            Z(:,1:(R(2)-1)) = [];
            Z = imrotate(Z,90);
        end
        SKEL = bwmorph(Z,'skeleton');
        
        
        
        
        
        skel = [];
        [skel(:,2) skel(:,1)] = find(SKEL);
        
        
      
        for p = 1:size(skel,1)
            sampX = squeeze(D(skel(p,1),skel(p,2),:,:,:));
            sampX = permute(sampX,[1 3 2]);
            sz = size(sampX);
            sampX = reshape(sampX,[prod(sz(1:2)) sz(3)]);
            sigX{cnt} = sampX;
           
            
            
            sampY = squeeze(DM(skel(p,1),skel(p,2),:,:,:));
            sampY = round(sampY);
            sigY{cnt} = categorical(sampY');
            
            
            cnt = cnt + 1;
        end
        
        
        
        imshow(Z,[]);
        drawnow
    end
        
end