function [pY] = runZoomSequence2(I,MAG_LOC_net,pY,boxSequence,zoomSequence,uY,sY,ITER)
    disp = false;
    try
        PSZ = 500;
        for k = 1:size(I,3)
            tmp = I(:,:,k);
            tmp = padarray(tmp,[PSZ PSZ],'replicate','both');
            newI(:,:,k) = tmp;
        end
        I = newI;
        pY = bsxfun(@plus,pY,[PSZ PSZ]);
        
        tm = clock;
        for p = 1:size(pY,1)
            cnt = 1;
            for z = 1:size(boxSequence,1)
                for iter = 1:ITER(z)
                    BOX = point2Box((double(round(pY(p,:,cnt)))),boxSequence(z,:));
                    [subI] = mimcrop(I,BOX,zoomSequence(z),[]);

                   



                    yPre = MAG_LOC_net.predict(subI(:,:,1));
                    %yPre = bsxfun(@times,yPre,sY);
                    %yPre = bsxfun(@plus,yPre,uY);
                    POS = flipdim(yPre(1:2),2);
                    if disp
                        imshow(subI,[]);
                        hold on
                        plot(POS(2),POS(1),'g*');
                        pause(.5)
                    end
                    %waitforbuttonpress
                    %POS = fitbitp(subI,NB,map{z},beta{z},WINDOW{z});
                    deltaZ = POS - [size(subI,2) size(subI,1)]/2;
                    pY(p,:,cnt+1) = pY(p,:,cnt) + zoomSequence(z)^-1*((deltaZ));
                    cnt = cnt + 1;
                end
            end
        end
        
        
        pY = bsxfun(@minus,pY,[PSZ PSZ]);
        
        etm = etime(clock,tm);
        fprintf(['Done running zoom sequence in:' num2str(etm) '\n']);
    catch ME
        getReport(ME)
    end
end