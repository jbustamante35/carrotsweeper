function [] = stripComponents1(fileName,obj,obj2nd,outName,backGround_oPath,qrcode_oPath,colorSample_oPath,extraData_oPath,conetainerData_oPath,disp)
            % read the image
            I = double(imread(fileName))/255;
            % split the image into bio and meta
            [dataStrip,bioStrip,cropLine,msg,qrCropBox] = splitMaizeSeedlingImage(I,20);
            % generate background signals
            [vB1,hB1,vB2,hB2,TOT,MSK] = generateBKsignals_wrapper(I,obj);
            % generate backgrouns mask
            [MSK] = generateBackgroundMask(I,TOT);
            
            [bioMask] = generateBackgroundMask(bioStrip,TOT((cropLine+1):end,:,:));
            [dataMask] = generateBackgroundMask(dataStrip,TOT(1:cropLine,:,:));
            
            
            dataMask(1,:) = 0;
            dataMask(end,:) = 0;
            dataMask = ~imfill(~dataMask,fliplr(round(qrCropBox(1:2) + 10)),8);
            dataMask = bwareaopen(dataMask,1000);
            dataMask = imclearborder(dataMask);
            extraR = regionprops(dataMask);
            for e = 1:numel(extraR)
                extraR(e).BoundingBox(1:2) = extraR(e).BoundingBox(1:2) - 101;
                extraR(e).BoundingBox(3:4) = extraR(e).BoundingBox(3:4) + 2*101;
                boundingBox = extraR(e).BoundingBox(3:4);
                extraOne = imcrop(dataStrip,extraR(e).BoundingBox);
                extraMask = imcrop(dataMask,extraR(e).BoundingBox);
                save([extraData_oPath outName '_' num2str(e) '.mat'],'extraOne','extraMask','boundingBox');
                %imshow(extraOne,[]);
                %waitforbuttonpress
                %here = 1;
            end
            
            
            
            agumentBioMask = zeros(size(dataMask,1),size(dataMask,2));
            agumentBioMask = [agumentBioMask;bioMask];
            bidx = find(agumentBioMask);
            colorSample = [];
            for k = 1:3
                tmp = I(:,:,k);
                colorSample = [colorSample tmp(bidx)];
            end
            
            SNIP = 100;
            if ~isempty(obj2nd)
                cidx = obj2nd.cluster(colorSample);
                Z = zeros(size(agumentBioMask,1),size(agumentBioMask,2));
                Z(bidx) = cidx;
                Z((end-SNIP):end,:) = 0;
                %containerMask = Z ~= 0 & Z ~= 4 & Z ~= 6 & Z ~= 2 & Z ~= 5;
                containerMask = (Z == 2 | Z == 3 | Z == 4 | Z == 5) & Z ~= 6;
                containerMask = bwareaopen(logical(containerMask),9000);
                containerMask = imfill(containerMask,'holes');
                continerR = regionprops(logical(containerMask),'MinorAxisLength','BoundingBox');
                continerR([continerR.MinorAxisLength] < 300) = [];
                coneBBshift = [71 101];
                for cc = 1:numel(continerR)
                    continerR(cc).BoundingBox(1:2) = continerR(cc).BoundingBox(1:2) - coneBBshift;
                    continerR(cc).BoundingBox(3) = continerR(cc).BoundingBox(3) + 2*coneBBshift(1);
                    continerR(cc).BoundingBox(4) = continerR(cc).BoundingBox(4) + coneBBshift(2);
                    
                    
                    coneImage = imcrop(I,continerR(cc).BoundingBox);
                    coneMask = imcrop(containerMask,continerR(cc).BoundingBox);
                    coneLabel = imcrop(Z,continerR(cc).BoundingBox);
                    boundingBox = continerR(cc).BoundingBox;
                    save([conetainerData_oPath outName '_' num2str(cc) '.mat'],'coneImage','coneMask','boundingBox','coneLabel');
                    %imshow(bsxfun(@times,coneImage,coneMask),[]);
                    %drawnow
                    %waitforbuttonpress
                end
                
                
            end
          
            
            qrCropBox(1:2) = qrCropBox(1:2) - 101;
            qrCropBox(3:4) = qrCropBox(3:4) + 2*101;
            QR = imcrop(I,qrCropBox);
            MQR = imcrop(MSK,qrCropBox);
            MQR = bwlarge(MQR);
           
            save([backGround_oPath 'v' filesep '1' filesep outName '.mat'],'vB1');
            save([backGround_oPath 'h' filesep '1' filesep outName '.mat'],'hB1');
            save([backGround_oPath 'v' filesep '2' filesep outName '.mat'],'vB2');
            save([backGround_oPath 'h' filesep '2' filesep outName '.mat'],'hB2');
            
            boundingBox = qrCropBox;
            save([qrcode_oPath outName '.mat'],'QR','MQR','boundingBox');
            
            save([colorSample_oPath outName '.mat'],'colorSample');
            
            
            if disp
                
                out = flattenMaskOverlay(I,~logical(MSK),.5,'r');

                imshow(out,[]);
                drawnow
                imshow(bsxfun(@times,QR,MQR),[]);
                drawnow
            end
end