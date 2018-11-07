function [out] = isolateRoots_ver2(fileName,disp)
    try
        SNIP = 50;
        
        
        I = double(imread(fileName))/255;
        I(:,end) = [];
        I = handleFLIP(I,[]);
        
        para.scales.value = 7;
        para.resize.value = 1;
        K = surKur(I,para);
        
        para.vars.sig.value = 2;
        para.vars.gradPara.sz.value = 2;
        para.vars.gradPara.method.value = 'finite';
        cM = cornerMap(I,para);
        
        cM = bindVec(cM);
        tipMask = cM > graythresh(cM);
        
        rM = bindVec(K(:,:,1));
        rootMask = rM > graythresh(rM);
        
        eM = bindVec(K(:,:,2));
        edgeMask = eM < graythresh(eM);
       
        fidx = find(eM == rM);
        rM(fidx) = 0;
        
        
        [g1 g2] = gradient(rM);
        
        
         
            t1 = ba_interp2(g1,c(2),r(2));
            t2 = ba_interp2(g2,c(2),r(2));
            
             
            T = [t2 t1];
            T = T / norm(T);
            N = [T(2) -T(1)];
            
            initD = [T;N];
            
            %quiver(curve(e).data(1,:),curve(e).data(2,:),d1X1,d1X2)
            %quiver(curve(e).data(1,:),curve(e).data(2,:),-d1X2,d1X1,'r')
            
            
        
        B = double(bwdist(rootMask));
            
            midlines = trackFromPointAtGradient(rM,[c(2);r(2)],initD,5000,20,pi/4,[20 200],.01);
                                
        
        
        
        
        
        
        [rd1 rd2] = find(rootMask);
        XD = [rd1';rd2'];
        rd3 = ba_interp2((rM),XD(2,:),XD(1,:));
        XD = [XD;rd3];
        T = Radjacency(XD,50);
        
        %[c r V] = impixel(I);
        [str] = snapTo([rd1 rd2],[r(1) c(1)]);
        [stp] = snapTo([rd1 rd2],[r(2) c(2)]);
        [p pc] = dijkstra(T,str,stp);
        gamma = XD(1:2,p);
        
        
        
        
        Io = I;
        Is = imresize(I,.25);
        BK = imclose(Is,strel('disk',21));
        BK = imfilter(BK,fspecial('disk',21),'replicate');
        BK = imresize(BK,size(Io));
        I = I - BK;
        I = I - min(I(:));
        I = I / max(I(:));
        I = imfilter(I,fspecial('gaussian',21,6),'replicate');
        [g1 g2] = gradient(I);
        thresh = graythresh(I);
        skel = I < thresh;
        skel = bwareaopen(skel,500);
        cB = imclearborder(skel);
        skel = skel & ~ cB;
        
        
        %%%%%%%%%%%%%%%
        %{
        Io = I;
        I = [I(:,1:SNIP) I];
        h = fspecial('gaussian',11,11);
        I = imfilter(I,h,'replicate');
        [g1 g2] = gradient(I);
        d = (g1.^2 + g2.^2).^.5;
        thresh = graythresh(d);
        MASK = d > thresh;
        skel = bwmorph(MASK,'skel',inf);
        skel = bwmorph(skel,'spur',inf);
        skel = skel(:,SNIP+1:end);
        skel(:,1) = 1;
        skel = imfill(skel,'holes');
        filler = sum(skel(:,1:2),2) == 2;
        skel(:,1) = filler;
        %}
        %%%%%%%%%%%%%%%






        
        C = contourc(double(skel),[1 1]);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
        % get the curve structure
        str = 1;
        c = 1;
        clear curve
        while str < size(C,2)
            ed = str + C(2,str);
            curve(c).level = C(1,str);
            curve(c).data = C(:,str+1:ed);
            curve(c).length = size(curve(c).data,2);
            c = c + 1;
            str = ed + 1;
        end

        ridx = find([curve.length] < 100);
        curve(ridx) = [];


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % find the highest curvature
        KSNIP = 50;
        SMOOTH_VALUE = 20;
        for e = 1:numel(curve)
            o = cwtK(curve(e).data',{SMOOTH_VALUE});
            [J tipIDX{e}] = min(o.K);
            o = cwtK(curve(e).data',{20});
            [J fine_tipIDX] = min((o.K(tipIDX{e}-KSNIP:tipIDX{e}+KSNIP)));
            tipIDX{e} = tipIDX{e} + (fine_tipIDX - KSNIP - 1);
        end

       
        
        B = double(bwdist(~skel));
        g1 = -g1;
        g2 = -g2;
        
        for e = 1:numel(curve)
            
            t1 = ba_interp2(g1,curve(e).data(1,tipIDX{e}),curve(e).data(2,tipIDX{e}));
            t2 = ba_interp2(g2,curve(e).data(1,tipIDX{e}),curve(e).data(2,tipIDX{e}));
            
             
            N = [t2 t1];
            N = -N / norm(N);
            T = [N(2) -N(1)];
            initD = [T;N];
            
            %quiver(curve(e).data(1,:),curve(e).data(2,:),d1X1,d1X2)
            %quiver(curve(e).data(1,:),curve(e).data(2,:),-d1X2,d1X1,'r')
            
            midlines(e).data = trackFromPointAtGradient(B,curve(e).data(:,tipIDX{e}),initD,5000,20,pi/2,[20 200],.3);
                                
        end

        if disp
            hold off
            imshow(Io,[]);
            hold on
            for e = 1:numel(curve)
                plot(curve(e).data(1,:),curve(e).data(2,:),'r');
                plot(midlines(e).data(1,:),midlines(e).data(2,:),'g');
            end

            for e = 1:numel(curve)
                plot(curve(e).data(1,tipIDX{e}),curve(e).data(2,tipIDX{e}),'g*');
            end
            drawnow
        end
        out.midlines = midlines;
        out.contours = curve;
    catch ME;
        out.ME = ME;
    end
end