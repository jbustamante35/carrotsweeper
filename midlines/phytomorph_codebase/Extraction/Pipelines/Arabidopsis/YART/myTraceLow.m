function [pheno] = myTraceLow(I,X,disp,toC)
    %[x1 x2 V] = impixel(I);    
    x1 = X(:,1);
    x2 = X(:,2);
    
    
    
    if toC 
        I = imcomplement(I);
    end
    % get surface kurva;ture
    para.scales.value=5;    
    para.resize.value=1;
    K = surKur(I,para);
    % threshold the stack on curvature
    tmp = bindVec(K(:,:,1));
    B = tmp > graythresh(tmp);
    SKEL = bwmorph(B,'skel',inf);
    for k = 1:5
        ep = bwmorph(SKEL,'endpoints');
        SKEL(find(ep)) = 0;
    end
   
    [x y] = find(SKEL);
    DP = [x y]';
    T = Radjacency(DP,3);
    
    
    
   
    x1 = reshape(x1,[2 size(x1,1)/2]);
    x2 = reshape(x2,[2 size(x2,1)/2]);
    for p = 1:size(x1,2)
        idx = [];
        for e = 1:size(x1,1)
            delta = (x - x2(e,p)).^2 + (y - x1(e,p)).^2;
            [~,idx(e)] = min(delta);
        end
        [path , pathcost]  = dijkstra(T , idx(1) , idx(2));
        pheno(p).gamma = DP(:,path);
    end
    
    
    if disp
        imshow(I,[])
        imshow(cat(3,SKEL,B,tmp),[]);
        hold on
        for e = 1:numel(pheno)
            plot(pheno(e).gamma(2,:),pheno(e).gamma(1,:),'r')
        end
    end


end

function [] = generateCropBox(I,X,EXTRA)
    UL = X(1,:) - EXTRA;
    HW = X(2,:) - X(1,:);
    HW = HW + 2*EXTRA;
end