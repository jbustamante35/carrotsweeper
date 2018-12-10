function [] = viewData(m,d,squareFlag,measureFlag,newCompute)
    mov = figure;
    if ~isempty(d)
        dat = figure;
    end
    if ischar(m)
        tmp = load(m);
        m = tmp.miniStack/255;
    end
    
    if squareFlag
        [m,mm] = squareMini(m,20);
        m = m * 255;
        m = bsxfun(@times,m,cat(3,mm,mm,mm));
    end
    
    if measureFlag
        [O] = measureMovie(m,mm,0);
    end
    
    if newCompute
        sz = size(m);
        for e = 1:size(m,4)
            m(:,:,:,e) = bsxfun(@times,m(:,:,:,e),mm);
            %hsv(:,:,:,e) = rgb2hsv_fast(m(:,:,:,e));
        end
        dm = imerode(mm,strel('disk',25,0));
        h = squeeze(hsv(:,:,1,:));
        h = reshape(h,[prod(sz(1:2)) sz(4)]);
        %m = reshape(m,[prod(sz(1:3)) sz(4)]);
        %mmvec = mm(:);
        sh = h(find(dm),:);
        srm
    end
    
    
    for e = 1:1:size(m,4)
        figure(mov);
        imshow(m(:,:,:,e),[]);
        
        if ~isempty(d)
            figure(dat);
            plot(d','k');
            hold on
            plot(e,d(1,e),'ro')
            hold off
            drawnow
        end
        
        title(num2str(e))
        drawnow
        if any(e == [55 75 85 105])
        %    waitforbuttonpress
        end
        %pause(.3);
    end
end