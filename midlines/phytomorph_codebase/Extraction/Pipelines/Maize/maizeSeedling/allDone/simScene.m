function [] = simScene(bkFileList,qrFileList,coneFileList,imgSZ,oSIM,R,disp)

    try
        toQR = logical((randi(2)-1));
        toCONE = logical((randi(2,[1 3])-1));

        n1 = randi(numel(bkFileList{1}));
        n2 = randi(numel(bkFileList{2}));
        n3 = randi(numel(bkFileList{3}));
        n4 = randi(numel(bkFileList{4}));

        q1 = randi(numel(qrFileList));


        c1 = randi(numel(coneFileList{1}));
        c2 = randi(numel(coneFileList{2}));
        c3 = randi(numel(coneFileList{3}));
        cName1 = coneFileList{1}{c1};
        cName2 = coneFileList{2}{c2};
        cName3 = coneFileList{3}{c3};
        
        cone1 = load(coneFileList{1}{c1},'coneImage','boundingBox','coneLabel','coneMask');
        cone2 = load(coneFileList{2}{c2},'coneImage','boundingBox','coneLabel','coneMask');
        cone3 = load(coneFileList{3}{c3},'coneImage','boundingBox','coneLabel','coneMask');
        cone1.coneMask = imclose(cone1.coneMask,strel('square',121));
        cone1.coneMask = imfill(cone1.coneMask,'holes');
        cone2.coneMask = imclose(cone2.coneMask,strel('square',121));
        cone2.coneMask = imfill(cone2.coneMask,'holes');
        cone3.coneMask = imclose(cone3.coneMask,strel('square',121));
        cone3.coneMask = imfill(cone3.coneMask,'holes');
        cone1.coneImage = bsxfun(@times,cone1.coneImage,cone1.coneMask);
        cone2.coneImage = bsxfun(@times,cone2.coneImage,cone2.coneMask);
        cone3.coneImage = bsxfun(@times,cone3.coneImage,cone3.coneMask);
        % set Y location for cone tainers

        Yloc = max([cone1.boundingBox(2)  cone2.boundingBox(2)  cone3.boundingBox(2)])+100;
        Yloc = imgSZ(1) - min([cone1.boundingBox(4) cone2.boundingBox(4) cone3.boundingBox(4)]);
        cone1.boundingBox(2) = Yloc;
        cone2.boundingBox(2) = Yloc;
        cone3.boundingBox(2) = Yloc;

        % set X location
        initX = 200;
        Xloc = linspace(initX,imgSZ(2)-10*initX,3);
        deltaX = randi(150,1,3);
        Xloc = Xloc + deltaX;
        cone1.boundingBox(1) = Xloc(1);
        cone2.boundingBox(1) = Xloc(2);
        cone3.boundingBox(1) = Xloc(3);



        QR = load(qrFileList{q1},'QR','boundingBox','MQR');
        QRLOC(1) = randi(imgSZ(2) - 500);
        QRLOC(2) = randi(round(randi(imgSZ(1)/2)));
        QR.boundingBox(1:2) = QRLOC;


        load(bkFileList{1}{n1},'hB1')
        load(bkFileList{2}{n2},'hB2')
        load(bkFileList{3}{n3},'vB1')
        load(bkFileList{4}{n4},'vB2')

        [simBK] = generateBackground(vB1,hB1,vB2,hB2,imgSZ);



        
        
        TOP_blend = zeros(size(simBK,1),size(simBK,2));
        TOP_color = zeros(size(simBK));

        countType = [0 0];
        if toCONE(1)
            cone1.toDropMask = cone1.coneLabel ~= 0 &  cone1.coneLabel ~= 2 &  cone1.coneLabel ~= 0;
            cone1.toDropMask = imfill(cone1.toDropMask,'holes');
            cone1.toDropMask = cone1.coneMask;
            if contains(cName1,'without')
                countType(1) = countType(1) + 1;
            else
                countType(2) = countType(2) + 1;
            end
            [TOP_color,TOP_blend] = dropImage(TOP_color,cone1.boundingBox,cone1.coneImage,cone1.toDropMask,TOP_blend);
        end


        if toCONE(2)
            cone2.toDropMask = cone2.coneLabel ~= 0 &  cone2.coneLabel ~= 2 &  cone2.coneLabel ~= 0;
            cone2.toDropMask = imfill(cone2.toDropMask,'holes');
            cone2.toDropMask = cone2.coneMask;
             if contains(cName2,'without')
                countType(1) = countType(1) + 1;
            else
                countType(2) = countType(2) + 1;
            end
            [TOP_color,TOP_blend] = dropImage(TOP_color,cone2.boundingBox,cone2.coneImage,cone2.toDropMask,TOP_blend);
        end

        if toCONE(3)
            cone3.toDropMask = cone3.coneLabel ~= 0 &  cone3.coneLabel ~= 2 &  cone3.coneLabel ~= 0;
            cone3.toDropMask = imfill(cone3.toDropMask,'holes');
            cone3.toDropMask = cone3.coneMask;
             if contains(cName3,'without')
                countType(1) = countType(1) + 1;
            else
                countType(2) = countType(2) + 1;
            end
            [TOP_color,TOP_blend] = dropImage(TOP_color,cone3.boundingBox,cone3.coneImage,cone3.toDropMask,TOP_blend);
        end

        if toQR
            [TOP_color,TOP_blend,QRLOC,wholeQRMask] = dropImage(TOP_color,QR.boundingBox,QR.QR,QR.MQR,TOP_blend);
        else
            wholeQRMask = NaN;
        end


        TOP_blend_blur = imfilter(TOP_blend,fspecial('disk',11),'replicate');
        TOP_blend = TOP_blend.*TOP_blend_blur;
        BOT_blend = 1 - TOP_blend;


        TOP_color = bsxfun(@times,TOP_color,TOP_blend);
        BOT_color = bsxfun(@times,simBK,BOT_blend);
        TOT_color = BOT_color + TOP_color;
        %imshow(TOT_color,[]);


        numCONE = sum(toCONE);
        imwrite(TOT_color,[oSIM num2str(R) '.tif'])
        save([oSIM num2str(R) '.mat'],'toQR','numCONE','wholeQRMask','TOP_blend','countType');

        if disp
            imshow(TOT_color,[])
            countType
            drawnow
        end
        
        
    catch ME
        ME
    end
    %imshow(simBK,[]);
    %drawnow