function [X,Y] = sampleSparkleJitterSequence(FileList,sparkleList,pointList,boxSequence,zoomSequence,jiggleRad,jiggleNumber,rotPara,disp)
   
    % if the same box - then replicate
    if size(boxSequence,1) == 1
        boxSequence = repmat(boxSequence,[numel(zoomSequence) 1]);
    end
    
    % theta vector
    TH = linspace(-pi,pi,100);
    
    cnt = ones(1,size(boxSequence,1));
    
    fprintf(['Starting image jitter sample for sparkle sequence\n']);
    % for each image
    for e = 1:numel(FileList)
        close all
        I = double(imread(FileList{e}))/255;
        
        if disp
            h1 = figure;
            subplot(1,2,1)
            imshow(I,[]);
            hold on
        end
        
        % for each point
        for p = 1:size(pointList{e},1)
        %p = e;
            % for each zoom amount
            for z = 1:size(boxSequence,1)
                % for some number of jiggles
                for j = 1:jiggleNumber
                    
                    randRad = randi(jiggleRad(z),1);
                    randTH = TH(randi(numel(TH),1));

                    delta = [randRad*cos(randTH) randRad*sin(randTH)];
                    %delta = 0;
                    rndP = pointList{e}(p,:) + delta;
                    
                    BOX = point2Box(rndP,boxSequence(z,:));
                    
                    
                    [subI,rotM,disU] = mimcrop(I,BOX,zoomSequence(z),rotPara);
                    %[subI,rotM,disU] = mimcrop(I,BOX,zoomSequence(z),[]);
                    
                    
                    tmpY = sparkleList{e};
                    %tmpY = [0 399] + rndP;
                    %tmpU = rndP;
                    %tmpU = mean(tmpY,1);
                    
                    
                    tmpY = bsxfun(@minus,tmpY,disU);
                    tmpY = (rotM'*tmpY')';
                    %tmpY = bsxfun(@plus,tmpY,tmpU);
                    
                    %tmpY = bsxfun(@plus,tmpY,disU);
                    
                    tmpY = tmpY * zoomSequence(z);
                    tmpY = bsxfun(@plus,tmpY,[size(subI,2) size(subI,1)]/2);
                    
                    
                    
                    if cnt(z) == 1
                         X{z} = zeros([size(subI) jiggleNumber*size(pointList{e},1)*numel(FileList)]);
                         Y{z} = zeros([2 size(tmpY,1) jiggleNumber*size(pointList{e},1)*numel(FileList)]);
                    end
                    
                    
                    X{z}(:,:,:,cnt(z)) = subI;
                    
                    Y{z}(:,:,cnt(z)) = tmpY';
                    
                    
                    cnt(z) = cnt(z) + 1;
                    
                    if disp
                        
                        figure(h1)
                        subplot(1,2,1)
                        rectangle('Position',BOX,'EdgeColor','g');
                        subplot(1,2,2)
                        imshow(subI,[]);
                        hold on
                        plot(Y{z}(1,:,cnt(z)-1),Y{z}(2,:,cnt(z)-1),'g.')
                        hold off
                        drawnow
                    end
                end


                

            end
            
            
        %end
    end
end