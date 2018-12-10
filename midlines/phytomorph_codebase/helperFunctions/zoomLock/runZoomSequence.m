function [pY] = runZoomSequence(I,NB,pY,boxSequence,zoomSequence,map,beta,WINDOW,ITER)
    tm = clock;
    for p = 1:size(pY,1)
        cnt = 1;
        for z = 1:numel(map)
            for iter = 1:ITER(z)
                BOX = point2Box((double(round(pY(p,:,cnt)))),boxSequence(z,:));
                [subI] = mimcrop(I,BOX,zoomSequence(z),[]);
                
                POS = fitbitp(subI,NB,map{z},beta{z},WINDOW{z});
                deltaZ = POS - [size(subI,2) size(subI,1)]/2;
                pY(p,:,cnt+1) = pY(p,:,cnt) + zoomSequence(z)^-1*((deltaZ));
                cnt = cnt + 1;
            end
        end
    end
    etm = etime(clock,tm);
    fprintf(['Done running zoom sequence in:' num2str(etm) '\n']);
end