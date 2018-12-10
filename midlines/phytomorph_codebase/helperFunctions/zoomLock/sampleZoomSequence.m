function [X,Y,PN] = sampleZoomSequence(FileList,pointList,boxSequence,zoomSequence,jiggleRad,jiggleNumber,rotPara,disp)
   
    % if the same box - then replicate
    if size(boxSequence,1) == 1
        boxSequence = repmat(boxSequence,[numel(zoomSequence) 1]);
    end
    
    % theta vector
    TH = linspace(-pi,pi,100);
    
    cnt = ones(1,size(boxSequence,1));
    
    
   
    
    % for each image
    for e = 1:numel(FileList)
        e
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

            % for each zoom amount
            for z = 1:size(boxSequence,1)

                % for some number of jiggles
                for j = 1:jiggleNumber


                    randRad = randi(jiggleRad(z),1);
                    randTH = TH(randi(numel(TH),1));

                    delta = [randRad*cos(randTH) randRad*sin(randTH)];

                    rndP = pointList{e}(p,:) + delta;
                    
                    BOX = point2Box(rndP,boxSequence(z,:));
                    
                    
                    [subI] = mimcrop(I,BOX,zoomSequence(z),rotPara);
                    
                    
                    
                    if cnt(z) == 1
                        
                        
                         PN{z} = [];
                         X{z} = zeros([size(subI) jiggleNumber*size(pointList{e},1)*numel(FileList)]);
                         Y{z} = zeros([2 jiggleNumber*size(pointList{e},1)*numel(FileList)]);
                    end
                    
                    
                    X{z}(:,:,:,cnt(z)) = subI;
                    Y{z}(:,cnt(z)) = .5*[size(subI,2) size(subI,1)] - delta*zoomSequence(z);
                    cnt(z) = cnt(z) + 1;
                    
                    PN{z} = [PN{z};p];
                    
                    
                    if disp
                        
                        figure(h1)
                        subplot(1,2,1)
                        rectangle('Position',BOX,'EdgeColor','g');
                        subplot(1,2,2)
                        imshow(subI,[]);
                        hold on
                        plot(Y{z}(1,cnt(z)-1),Y{z}(2,cnt(z)-1),'g*')
                        hold off
                        drawnow
                    end
                end


                

            end
            
            
        end
    end
end