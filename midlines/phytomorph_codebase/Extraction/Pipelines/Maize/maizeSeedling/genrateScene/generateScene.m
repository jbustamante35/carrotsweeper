function [newImage] = generateScene(masterTable,qrDistribution,conetainerDistribution,sceneSize)


    simQRcode = 1;
    simConetainer = 1;
    
    numContainer = randi(2,1,3)-1;
    
    
    
    
    hBK = find(strcmp(masterTable.type,'background_h'));
    vBK = find(strcmp(masterTable.type,'background_v'));
    selH = randi(numel(hBK),1);
    selV = randi(numel(vBK),1);
    hBKimage = double(imread(masterTable.imageLocation{hBK(selH)}))/255;
    vBKimage = double(imread(masterTable.imageLocation{vBK(selV)}))/255;
    hBKimage = imresize(hBKimage,[1 sceneSize(2)]);
    hBKimage = repmat(hBKimage,[sceneSize(1) 1 1]);
    vBKimage = imresize(vBKimage,[sceneSize(1) 1]);
    vBKimage = repmat(vBKimage,[1 sceneSize(2) 1]);
    simBackground = hBKimage + vBKimage;
    
    
    if simQRcode
        qrFile = find(strcmp(masterTable.type,'qrObject'));
        selQR = randi(numel(qrFile),1);
        qrImage = double(imread(masterTable.imageLocation{qrFile(selQR)}))/255;
        
        qrLocation = (round(qrDistribution.random));
        %qrLocation = [masterTable.boundingBox1(qrFile(selQR)) masterTable.boundingBox2(qrFile(selQR))];
        
        newImage = blendEdge(simBackground,qrImage,qrLocation,21,21,[51 51]);
        
    end
    
    
    %numContainer(3) = 1;
    if simConetainer
        for e = 1:numel(numContainer)

            if numContainer(e)
            

                conetainerFile = find(strcmp(masterTable.type,'conetainer_whole_masked'));
                selCone = randi(numel(conetainerFile),1);
                rowIndex = conetainerFile(selCone);
                containerImage = double(imread(masterTable.imageLocation{rowIndex}))/255;



                IDX = 0;
                while IDX ~= e
                    [containerLocation IDX] = conetainerDistribution.random;
                end
                
                containerLocation = round((containerLocation));

                
                

                %containerLocation = [masterTable.boundingBox2(rowIndex) masterTable.boundingBox1(rowIndex)];
                containerLocation(1) = sceneSize(1) - size(containerImage,1);



                newImage = blendSimular(newImage,containerImage,containerLocation,11,11,[31 31]);
                %newImage = blendEdge(newImage,containerImage,containerLocation,5,3,[11 11]);
                %imshow(newImage,[]);
                %waitforbuttonpress
                
            end
        end
        
        
    end
    
end

function [newImage] = blendEdge(bkImage,objectImage,objectLocation,erodeA,filterA,filterSZ)
    try
        newImage = bkImage;
        blendImage = zeros(size(bkImage));

        % if object is too tall
        if objectLocation(1) + size(objectImage,1) > size(bkImage,1)
            toRemove = size(bkImage,1) - objectLocation(1) + size(objectImage,1);
            objectImage(end-toRemove,:,:) = [];
        end


        % if object is too wide
        if objectLocation(2) + size(objectImage,2) > size(bkImage,2)
            toRemove = size(bkImage,2) - objectLocation(2) + size(objectImage,2);
            objectImage(:,end-toRemove,:) = [];
        end



        blendImage(objectLocation(1):objectLocation(1)+size(objectImage,1)-1,...
                   objectLocation(2):objectLocation(2)+size(objectImage,2)-1,:) = 1;




        blendImage = imerode(blendImage,strel('disk',erodeA,0));
        blendImage = imfilter(blendImage,fspecial('gaussian',filterSZ,filterA),'replicate');


        newImage(objectLocation(1):objectLocation(1)+size(objectImage,1)-1,...
                   objectLocation(2):objectLocation(2)+size(objectImage,2)-1,:) = objectImage;

        newImage = blendImage.*newImage + (1-blendImage).*bkImage;
    catch ME
        ME
    end

end


function [newImage] = blendSimular(bkImage,objectImage,objectLocation,erodeA,filterA,filterSZ)
    try

        newImage = bkImage;
        blendImage = zeros(size(bkImage));

        % if object is too tall
        if objectLocation(1) + size(objectImage,1) > size(bkImage,1)
            toRemove = size(bkImage,1) - objectLocation(1) + size(objectImage,1);
            objectImage((end-toRemove):end,:,:) = [];
        end


        % if object is too wide
        if objectLocation(2) + size(objectImage,2) > size(bkImage,2)
            toRemove = size(bkImage,2) - objectLocation(2) + size(objectImage,2);
            objectImage(:,(end-toRemove):end,:) = [];
        end



        newImage(objectLocation(1):objectLocation(1)+size(objectImage,1)-1,...
                   objectLocation(2):objectLocation(2)+size(objectImage,2)-1,:) = objectImage;
        ridx = find(newImage==0);
        newImage(ridx) = bkImage(ridx);


        %{
        deltaI = sum((newImage - bkImage).^2,3).^.5;

        deltaPatch = deltaI(objectLocation(1):objectLocation(1)+size(objectImage,1)-1,...
                   objectLocation(2):objectLocation(2)+size(objectImage,2)-1,:);
        blendMask = zeros(size(blendImage));
        blendMask(objectLocation(1):objectLocation(1)+size(objectImage,1)-1,...
                   objectLocation(2):objectLocation(2)+size(objectImage,2)-1,:) = 1;

        blendImage = deltaI > graythresh(deltaPatch);

        blendImage = blendImage.*blendMask;
        blendSim = zeros(size(blendImage));
        blendSim(find(blendImage==0)) = .05;
        blendSim = blendSim.*blendMask;
        blendSim = imfilter(blendSim,fspecial('gaussian',filterSZ,filterA),'replicate');
        blendImage(find(blendImage==0)) = blendSim(find(blendImage==0));




        blendImage = blendImage.*blendMask;




        newImage = blendImage.*newImage + (1-blendImage).*bkImage;
        %}
    catch ME
        ME
    end
end

 
    
    
    
    