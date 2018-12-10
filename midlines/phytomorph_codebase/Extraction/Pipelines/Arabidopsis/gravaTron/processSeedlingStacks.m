function [SHOOT,ROOT] = processSeedlingStacks(cI)
r = [];
s = [];
close all
OFFSET = 100;
for e = 1:numel(cI)
    SHOOT(e).X = [];
    SHOOT(e).Y = [];
    ROOT(e).X = [];
    ROOT(e).Y = [];
    [bk dM] = getBKandDM(cI{e},20);
    [G] = makeMovieGray(cI{e});
    [d1 d2 d3] = gradient(G);
    for e = 1:size(G,3)
        
        para.scales.value = [10];
        para.resize.value = 1;
        [K] = surKur(G(:,:,e),para);
        K = K(:,:,1).*K(:,:,2).^-1;
        %K = imcomplement(bindVec(K));
       
        MSK = K > graythresh(K);
        MSK = imclearborder(MSK);
        imshow(MSK,[]);
        drawnow
    end
    
    %{
    sig = squeeze(mean(abs(dM).*abs(d3),2));
    ep = size(sig,2);
    for t = 1:size(sig,2)
        ts = abs(sig(:,t));
        ts = imfilter(ts,fspecial('disk',13));
        
        if t >= 5
            ep = s(e,t-1)+20;
        end
        
        
        
        [~,s(e,t)] = max(ts(1:ep));
        
        if t >= 5
            sp = r(e,t-1);
        else
            sp = s(e,t) + OFFSET;
            
        end
        sp = min(sp,numel(ts)-100);
        
        
        
        [~,r(e,t)] = max(gradient(ts(sp:end)));
        r(e,t) = r(e,t) + sp;
        
        
        
        imshow(cI{e}(:,:,:,t)/255,[])
        hold on
        SHOOT(e).X = [SHOOT(e).X;linspace(1,size(cI{e},2),2)];
        SHOOT(e).Y = [SHOOT(e).Y;s(e,t)*ones(1,2)];
        ROOT(e).X = [ROOT(e).X;linspace(1,size(cI{e},2),2)];
        ROOT(e).Y = [ROOT(e).Y;r(e,t)*ones(1,2)];
        
        plot(SHOOT(e).X(t,:),SHOOT(e).Y(t,:),'g')
        plot(ROOT(e).X(t,:),ROOT(e).Y(t,:),'r')
        hold off
        drawnow
    end
    %}
end
end