function [SP CA_map TOTMAP WF_C decay stopMap conc DTM CMF] = mSpeed(CCf,sz,tm_threshold,J,sJ,wv1,wv2,wv3)
    CA_map = zeros(sz(1),sz(2));
    TOTMAP = zeros(size(CA_map));
    stopMap = zeros(size(CA_map));
    SP = [];
    cnt = 1;
    decay = [];
    conc = [];
    DTM = [];
    CMF = {};
    for k = 1:CCf.NumObjects
        [w1 w2 w3] = ind2sub(sz,CCf.PixelIdxList{k});
        DT = max(w3) - min(w3);
        DT
        mean(w3)
        if all(w3 > tm_threshold) & numel(w1) > 200 & DT > 3
            UQ = unique(w3);
            DUR = 2*std(unique(w3));
            initf = UQ(1);
            
            
            
            for u = 1:numel(UQ)
                fidx = find(w3==UQ(u));
                GTH = [];
                for f = 1:numel(fidx)
                    TOTMAP(w1(fidx(f)),w2(fidx(f))) = 1;
                    try
                        GTH = [GTH;squeeze(sJ(w1(fidx(f)),w2(fidx(f)),UQ(u):(UQ(u)+149)))'];
                    catch
                        
                    end
                    %GTH = [GTH;[mean(wv1(w1(fidx(f)),w2(fidx(f)),UQ(u))) mean(wv2(w1(fidx(f)),w2(fidx(f)),UQ(u))) mean(wv3(w1(fidx(f)),w2(fidx(f)),UQ(u)))]];
                end
                if k == 1
                    CMF{u} = {};
                end
                try
                    CMF{u};
                catch
                    CMF{u} = [];
                end
                CMF{u} = [CMF{u};GTH];
            end
            
            UQ(1:2) = [];
            
            for u = 1:numel(UQ)
                Z = zeros(sz(1),sz(2));
                fidx = find(w3==UQ(u));
               
                for f = 1:numel(fidx)
                    Z(w1(fidx(f)),w2(fidx(f))) = 1;
                    
                    if u == numel(UQ)
                        stopMap(w1(fidx(f)),w2(fidx(f))) = 1;
                    end
                end
              
                if UQ(u) == 276
                    stop = 1;
                end
               
                
                R = regionprops(logical(Z),'Image','PixelIdxList','BoundingBox','Area','PixelIdxList');
                VS = [];
                for r = 1:numel(R)
                    if R(r).Area > 30
                        DIST = bwdist(~logical(R(r).Image));
                        wf = bwmorph(logical(R(r).Image),'thin',inf);
                        wfidx = find(wf);
                        v = max(DIST(:));
                        v = mean(DIST(wfidx));
                        [ww1 ww2] = find(wf);
                        WF_C{cnt} = bsxfun(@plus,[ww1 ww2],fliplr(R(r).BoundingBox(1:2)-.5));
                        cnt = cnt + 1;
                        if v == 4
                            stop = 0;
                        end
                        if v ~= inf
                            CA_map(R(r).PixelIdxList) = 2*v;
                        end
                        decay = [decay;[UQ(u)-initf 2*v]];
                        conc = [conc;[UQ(u)-initf mean(J(R(r).PixelIdxList))]];
                        
                        VS = [VS 2*v];
                        
                    end
                end
                
            end
            if any([R(r).Area] > 30)
                DTM = [DTM ;[initf mean(w3) DUR mean(VS)*DUR]];
            end
        end
    end
    fidx = find(CA_map~=0);
    SP = unique(CA_map(fidx));
    SP = CA_map(fidx);
end