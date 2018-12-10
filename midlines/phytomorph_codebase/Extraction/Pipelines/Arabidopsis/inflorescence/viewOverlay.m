function [] = viewOverlay(fileList,idxL,startPoint,initIDX,PATHS)
    for e = 1:numel(fileList)
        I = imread(fileList{e});
        MSK = zeros(size(I));
        MSK = bwareaopen(MSK,300);
        
        
        
        
            imshow(I,[]);
            hold on
            %{
        for sp = 1:size(startPoint,1)
        
            
            
            
        
            if e == 1
                curve = PATHS{sp}{e}{initIDX(sp)};
                plot(curve(:,2),curve(:,1),'m')
                
                
            end
        end
        %}
        
        
        
        fidx = idxL(:,2) == e;
        MSK(idxL(fidx,1)) = 1;
        SKEL = bwmorph(MSK,'skeleton',inf);
        [skel(:,1) skel(:,2)] = find(SKEL);
        for p = 1:size(startPoint,1)
            midx = snapTo(skel,startPoint(p,:));
        end
        
        
        EP = bwmorph(SKEL,'endpoints');
        ep = [];
        [ep(:,1) ep(:,2)] = find(EP);
        
        
         
        BP = bwmorph(SKEL,'branchpoints');
        bp = [];
        [bp(:,1) bp(:,2)] = find(BP);
        
        TP = [ep;bp];
        T = Radjacency(TP,2^.5+eps);
        
        
        imshow(I,[]);
        hold on
        plot(skel(:,2),skel(:,1),'b.','MarkerSize',2);
        plot(ep(:,2),ep(:,1),'go')
        plot(bp(:,2),bp(:,1),'ko')
        hold off
        
        
    end
end