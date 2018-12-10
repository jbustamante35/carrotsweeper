function [data] = graviFold(M,sz)
        sz = [0 sz 0];
        fidx = find(sz==0);
        for e = 1:numel(fidx)-1
           SZ{e} = sz(fidx(e)+1:fidx(e+1)-1); 
        end
        str = 1;
        stp = prod(SZ{1}(1:3));
        tmp1 = M(str:stp);
        str = stp + 1;
        stp = stp + prod(SZ{2}(1:3));
        tmp2 = M(str:stp);
        str = stp + 1;
        stp = stp + prod(SZ{3}(1:3));
        tmp3 = M(str:stp);
        str = stp + 1;
        stp = stp + prod(SZ{4}(1:3));
        tmp4 = M(str:stp);
        %{
        str = stp + 1;
        stp = stp + prod(SZ{5}(1:2));
        tmp5 = M(str:stp);
        str = stp + 1;
        stp = stp + prod(SZ{6}(1:2));
        tmp6 = M(str:stp);
        str = stp + 1;
        stp = stp + prod(SZ{7}(1:2));
        tmp7 = M(str:stp);
        %}
        tmp1 = reshape(tmp1,SZ{1}(1:3));
        tmp2 = reshape(tmp2,SZ{2}(1:3));
        tmp3 = reshape(tmp3,SZ{3}(1:3));
        tmp4 = reshape(tmp4,SZ{4}(1:3));
        
        for tm = 1:size(tmp1,3)
            MID{tm} = arcLength(tmp1(:,:,tm),'arcLen',0);
        end
        
        tmp6 = measureTipAngle(MID,10);
        tmp5 = measureLength(MID);
        tmp7 = measureKurvature(MID,30);
        %tmp7 = reshape(tmp7,SZ{7}(1:2));
        data.kernelContour = tmp2;
        data.upperRoot = tmp3;
        data.lowerRoot = tmp4;
        data.midline = tmp1;
        data.length = tmp5;
        data.angle = tmp6;
        data.kurvature = tmp7;
        
end



function [ta] = measureTipAngle(midline,N)
    for e =1:numel(midline)
        sz = size(midline{e},1);
        A = mean(diff(midline{e}(1:min(N,sz),:),1,1),1);
        ta(e) = -atan2(A(1),A(2));
    end
end


function [l] = measureLength(midline)
    for e = 1:numel(midline)
        l(e) = size(midline{e},1);
    end
end


function [k] = measureKurvature(midline,N)
    para{1} = 3;
    k = zeros(N,numel(midline));
    for e = 1:numel(midline)
        try
            [out] = cwtK_imfilter(midline{e}(1:N,:),para);
            k(:,e) = out.K;
        catch
            size(midline{e},1)
        end
    end
end