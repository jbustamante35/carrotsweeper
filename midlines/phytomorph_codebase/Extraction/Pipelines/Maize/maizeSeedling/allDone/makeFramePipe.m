function [BL] = makeFramePipe()
    close all
    W = 17;
    w = 4;
    B = zeros(2*W+1);
    B((end-1)/2+1,(end-1)/2+1) = 1;
    B((end-1)/2-w+1:(end-1)/2+w+1,1:((end-1)/2+1+w)) = 1;
    imshow(B,[]);

    BB = [];
    BB(:,:,1) = B;
    for e = 2:4
        BB(:,:,e) = imrotate(BB(:,:,e-1),90);
        imshow(BB(:,:,e),[]);
        drawnow
    end
    BL = [];
    cnt = 1;
    for e = 2:4
        v = nchoosek(1:4,e);
       v
        for k = 1:size(v,1)
            newB = zeros(size(B));
            for i = 1:size(v,2)
                newB = newB | (BB(:,:,v(k,i)) == 1);
            end
            newB = newB .* sum(newB(:)).^-1 + e*.00001;
            BL(:,:,cnt) = newB;

            imshow(BL(:,:,cnt),[])
            drawnow
            %waitforbuttonpress
            cnt = cnt + 1;

        end

    end
end