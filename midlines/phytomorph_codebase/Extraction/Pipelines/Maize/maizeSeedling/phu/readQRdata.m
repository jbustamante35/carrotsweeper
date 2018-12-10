function [msg] = readQRdata(I,points,offset,disp)
    try
        for e = 1:size(points,3)




            W(1) = points(2,1,e) - points(1,1,e);
            W(2) = points(5,1,e) - points(4,1,e);
            H(1) = points(4,2,e) - points(1,2,e);
            H(2) = points(5,2,e) - points(2,2,e);


            BOX = [points(1,:,e) + offset [mean(W) mean(H)]-2*offset];
            subQR = imcrop(I,BOX);



            if ~isempty(subQR)
                if size(subQR,1) > 1 & size(subQR,2) > 1
                    subQR = rgb2gray(subQR);
                    subQR = imadjust(subQR);
                    subQR = imsharpen(subQR,'Amount',2);
                    msg{e} = decode_qr(subQR);
                end
            end

            if disp
                imshow(subQR,[]);
                drawnow
                %waitforbuttonpress
            end

        end
    catch ME
        ME
    end
end