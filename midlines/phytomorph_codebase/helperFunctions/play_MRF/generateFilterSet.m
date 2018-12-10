function [baseFilter f1 f2] = generateFilterSet(T,P)
    pidx = find(P);
    pidx = 1:numel(P);
    ptms = permn(pidx,3);
    bidx = reshape(1:81,[9 9]);
    baseFilter = zeros(2^9,numel(pidx));
    f1 = zeros(2^9,numel(pidx),2^9);
    f2 = f1;
    for b = 1:numel(pidx)
        tmp = reshape(T(pidx(b),:),[3 3]);
        tmp(5) = 0;
        fidx = [5;find(tmp)];
        for e = 1:numel(fidx)
            z = zeros(3);
            z(fidx(e)) = 1;
            z = imresize(z,[9 9],'nearest');
            zidx = find(z);
            idx{e} = bidx(zidx);
        end
        
        for p = 1:size(ptms)
            z = zeros(9);
            for e = 1:size(ptms,2)
                tmp_mini = T(ptms(p,e),:);
                z(idx{e}) = tmp_mini;
            end
            %imshow(imresize(z,10,'nearest'),[]);
            %drawnow
            %swaitforbuttonpress
            val = vTest(z,P);
            baseFilter(ptms(p,1),b) = logical(baseFilter(ptms(p,1),b)) | logical(val);
            f1(ptms(p,2),b,ptms(p,1)) = logical(f1(ptms(p,2),b,ptms(p,1))) | logical(val);
            f2(ptms(p,3),b,ptms(p,1)) = logical(f2(ptms(p,3),b,ptms(p,1))) | logical(val);
            %{
            if vTest(z,P)
                z
                p
                ptms(p,1)
            end
            %}
        end
    end
end