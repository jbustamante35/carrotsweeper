function [OUT_SQUARE,IN_SQUARE,RED_SQUARE,RIGHT_STRIP_SQUARE,HEADER_SQUARE,CARROT_SQUARE,ORIN] = findDarkLines(I,n,oPath)



    angleSweep = 20;
    smallThresh = 100000;
    MinLength_blackFrame = 1000;
    MinLength_redFrame = 200;
    labI = rgb2lab(double(I)/255);
   
    
    
    darkMask = labI(:,:,1)/100 < graythresh(labI(:,:,1)/100);
    SZ = size(darkMask);
    darkMask = logical(imresize(darkMask,1/5));
    darkMask = bwareaopen(darkMask,round(smallThresh/5));
    darkMask = logical(imresize(darkMask,SZ));
    
    R_blackFrame = regionprops(logical(darkMask),'PixelIdxList');
    numberLines = numel(R_blackFrame)*4;
    
    redMask = labI(:,:,2) > 32;
    redMask = bwlarge(redMask,numel(R_blackFrame));
    R_redFrame = regionprops(logical(redMask),'PixelIdxList');
    
    
    blueMask = labI(:,:,3) < -10;
    blueMask = bwareaopen(blueMask,smallThresh/10);
    R_blueFrame = regionprops(logical(blueMask),'PixelIdxList');
    
    
    
    
    
    
    
    for e = 1:numel(R_blackFrame)
    
       
        tmpMask = zeros(size(darkMask));
        tmpMask(R_blackFrame(e).PixelIdxList) = 1;
        tmpMask = imfill(tmpMask,'holes');
        E = edge(tmpMask);
        dE = imdilate(E,strel('disk',15,0));
        
        
        
        
        % horizontal lines
        [H,T,RHO] = hough(E','Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[101 31]);
        Hlines{e} = houghlines(dE',T,RHO,P,'FillGap',10000,'MinLength',MinLength_blackFrame);


        % vectical lines
        [H,T,RHO] = hough(E,'Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[101 31]);
        Vlines{e} = houghlines(dE,T,RHO,P,'FillGap',10000,'MinLength',MinLength_blackFrame);
        

        %{
        imshow(E,[]);
        hold on
        for l = 1:numel(Hlines{e})
            L = [Hlines{e}(l).point1;Hlines{e}(l).point2];
            plot(L(:,2),L(:,1),'g');
            plot(Hlines{e}(l).point1(2),Hlines{e}(l).point1(1),'mo')
            plot(Hlines{e}(l).point2(2),Hlines{e}(l).point2(1),'mo')
            
            L = [Vlines{e}(l).point1;Vlines{e}(l).point2];
            plot(L(:,1),L(:,2),'g');
            plot(Vlines{e}(l).point1(1),Vlines{e}(l).point1(2),'mo')
            plot(Vlines{e}(l).point2(1),Vlines{e}(l).point2(2),'mo')
        end
        %}
        
        
        tmpMask = zeros(size(darkMask));
        tmpMask(R_blackFrame(e).PixelIdxList) = 1;
        tmpMask = imfill(tmpMask,'holes') - tmpMask;
        tmpMask = bwareaopen(tmpMask,smallThresh);
        
        SZ = size(tmpMask);
        tmpMask = logical(imresize(tmpMask,1/5));
        tmpMask = imclose(tmpMask,strel('disk',round(400/5)));
        tmpMask = logical(imresize(tmpMask,SZ));
        
        tmpMask = bwlarge(logical(tmpMask));
        E = edge(tmpMask);
        dE = imdilate(E,strel('disk',15,0));
        
        
        % horizontal lines
        [H,T,RHO] = hough(E','Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[201 31]);
        inner_Hlines{e} = houghlines(dE',T,RHO,P,'FillGap',10000,'MinLength',MinLength_redFrame);


        % vectical lines
        [H,T,RHO] = hough(E,'Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[201 31]);
        inner_Vlines{e} = houghlines(dE,T,RHO,P,'FillGap',10000,'MinLength',MinLength_redFrame);
        
        
        %%%%% BLUE
        tmpBlue = blueMask.*tmpMask;
        
        tmpBlue = logical(imresize(tmpBlue,1/5));
        tmpBlue = imclose(tmpBlue,strel('disk',round(101/5)));
        tmpBlue = logical(imresize(tmpBlue,size(tmpMask)));
        
        
        E = edge(tmpBlue);
        dE = imdilate(E,strel('disk',15,0));
        % horizontal lines
        [H,T,RHO] = hough(E,'Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[101 31]);
        blue_Vlines{e} = houghlines(dE,T,RHO,P,'FillGap',10000,'MinLength',MinLength_redFrame);
        
        
    end
    
    
    
    
    for e = 1:numel(R_redFrame)
    
       
        tmpMask = zeros(size(darkMask));
        tmpMask(R_redFrame(e).PixelIdxList) = 1;
        tmpMask = imfill(tmpMask,'holes');
        E = edge(tmpMask);
        dE = imdilate(E,strel('disk',11,0));
        
        
        
        
        % horizontal lines
        [H,T,RHO] = hough(E','Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[101 31]);
        Hlines_red{e} = houghlines(dE',T,RHO,P,'FillGap',10000,'MinLength',MinLength_redFrame);


        % vectical lines
        [H,T,RHO] = hough(E,'Theta',linspace(-angleSweep,angleSweep,50));
        P  = houghpeaks(H,2,'threshold',ceil(0.0001*max(H(:))),'NHoodSize',[101 31]);
        Vlines_red{e} = houghlines(dE,T,RHO,P,'FillGap',10000,'MinLength',MinLength_redFrame);
        

        %{
        imshow(E,[]);
        hold on
        for l = 1:numel(Hlines_red{e})
            L = [Hlines_red{e}(l).point1;Hlines_red{e}(l).point2];
            plot(L(:,2),L(:,1),'g');
            plot(Hlines_red{e}(l).point1(2),Hlines_red{e}(l).point1(1),'mo')
            plot(Hlines_red{e}(l).point2(2),Hlines_red{e}(l).point2(1),'mo')
            
            L = [Vlines_red{e}(l).point1;Vlines_red{e}(l).point2];
            plot(L(:,1),L(:,2),'g');
            plot(Vlines_red{e}(l).point1(1),Vlines_red{e}(l).point1(2),'mo')
            plot(Vlines_red{e}(l).point2(1),Vlines_red{e}(l).point2(2),'mo')
            
        end
        drawnow
        hold off
        %}
    end
    
    
    hL = [];
    for e = 1:numel(Hlines)
        HL = cat(3,flip([Hlines{e}(1).point1;Hlines{e}(1).point2],2),flip([Hlines{e}(2).point1;Hlines{e}(2).point2],2));
        VL = cat(3,([Vlines{e}(1).point1;Vlines{e}(1).point2]),([Vlines{e}(2).point1;Vlines{e}(2).point2]));
        L = cat(4,HL,VL);
        OUT_SQUARE{e} = squareStitch(I,L);
        
        
        HL = cat(3,flip([inner_Hlines{e}(1).point1;inner_Hlines{e}(1).point2],2),flip([inner_Hlines{e}(2).point1;inner_Hlines{e}(2).point2],2));
        VL = cat(3,([inner_Vlines{e}(1).point1;inner_Vlines{e}(1).point2]),([inner_Vlines{e}(2).point1;inner_Vlines{e}(2).point2]));
        L = cat(4,HL,VL);
        IN_SQUARE{e} = squareStitch(I,L);
        
        
        HL = cat(3,flip([Hlines_red{e}(1).point1;Hlines_red{e}(1).point2],2),flip([Hlines_red{e}(2).point1;Hlines_red{e}(2).point2],2));
        VL = cat(3,([Vlines_red{e}(1).point1;Vlines_red{e}(1).point2]),([Vlines_red{e}(2).point1;Vlines_red{e}(2).point2]));
        L = cat(4,HL,VL);
        RED_SQUARE{e} = squareStitch(I,L);
        
        
        VL = cat(3,([blue_Vlines{e}(1).point1;blue_Vlines{e}(1).point2]),([blue_Vlines{e}(2).point1;blue_Vlines{e}(2).point2]));
        HL = cat(3,flip([inner_Hlines{e}(1).point1;inner_Hlines{e}(1).point2],2),flip([inner_Hlines{e}(2).point1;inner_Hlines{e}(2).point2],2));
        L = cat(4,HL,VL);
        BLUE_STRIP_SQUARE{e} = squareStitch(I,L);
        
        VL = cat(3,BLUE_STRIP_SQUARE{e}(4:5,:),IN_SQUARE{e}(4:5,:));
        HL = cat(3,flip([inner_Hlines{e}(1).point1;inner_Hlines{e}(1).point2],2),flip([inner_Hlines{e}(2).point1;inner_Hlines{e}(2).point2],2));
        L = cat(4,HL,VL);
        LEFT_STRIP_SQUARE{e} = squareStitch(I,L);
        
        VL = cat(3,BLUE_STRIP_SQUARE{e}(2:3,:),IN_SQUARE{e}(2:3,:));
        HL = cat(3,flip([inner_Hlines{e}(1).point1;inner_Hlines{e}(1).point2],2),flip([inner_Hlines{e}(2).point1;inner_Hlines{e}(2).point2],2));
        L = cat(4,HL,VL);
        RIGHT_STRIP_SQUARE{e} = squareStitch(I,L);
        
        
        if all(inpolygon(RED_SQUARE{e}(:,1),RED_SQUARE{e}(:,2),LEFT_STRIP_SQUARE{e}(:,1),LEFT_STRIP_SQUARE{e}(:,2)))
            HEADER_SQUARE{e} = LEFT_STRIP_SQUARE{e};
            CARROT_SQUARE{e} = RIGHT_STRIP_SQUARE{e};
            ORIN(e) = true;
        else
            HEADER_SQUARE{e} = RIGHT_STRIP_SQUARE{e};
            CARROT_SQUARE{e} = LEFT_STRIP_SQUARE{e};
            ORIN(e) = false;
        end
        
    end
    
    
    hold off
    imshow(I,[]);
    hold on
    for e = 1:numel(Hlines)
        for l = 1:numel(Hlines{e})
            L = [Hlines{e}(l).point1;Hlines{e}(l).point2];
            plot(L(:,2),L(:,1),'k','LineWidth',3);
            plot(Hlines{e}(l).point1(2),Hlines{e}(l).point1(1),'k*')
            plot(Hlines{e}(l).point2(2),Hlines{e}(l).point2(1),'k*')

            L = [Vlines{e}(l).point1;Vlines{e}(l).point2];
            plot(L(:,1),L(:,2),'k','LineWidth',3);
            plot(Vlines{e}(l).point1(1),Vlines{e}(l).point1(2),'k*')
            plot(Vlines{e}(l).point2(1),Vlines{e}(l).point2(2),'k*')
        end
    end
    
    for e = 1:numel(inner_Hlines)
        for l = 1:numel(Hlines{e})
            L = [inner_Hlines{e}(l).point1;inner_Hlines{e}(l).point2];
            plot(L(:,2),L(:,1),'g','LineWidth',3);
            plot(inner_Hlines{e}(l).point1(2),inner_Hlines{e}(l).point1(1),'g*')
            plot(inner_Hlines{e}(l).point2(2),inner_Hlines{e}(l).point2(1),'g*')

            L = [inner_Vlines{e}(l).point1;inner_Vlines{e}(l).point2];
            plot(L(:,1),L(:,2),'g','LineWidth',3);
            plot(inner_Vlines{e}(l).point1(1),inner_Vlines{e}(l).point1(2),'g*')
            plot(inner_Vlines{e}(l).point2(1),inner_Vlines{e}(l).point2(2),'g*')
        end
    end
    
    for e = 1:numel(Hlines_red)
        for l = 1:numel(Hlines_red{e})
            L = [Hlines_red{e}(l).point1;Hlines_red{e}(l).point2];
            plot(L(:,2),L(:,1),'r','LineWidth',3);
            plot(Hlines_red{e}(l).point1(2),Hlines_red{e}(l).point1(1),'ro')
            plot(Hlines_red{e}(l).point2(2),Hlines_red{e}(l).point2(1),'ro')

            L = [Vlines_red{e}(l).point1;Vlines_red{e}(l).point2];
            plot(L(:,1),L(:,2),'r','LineWidth',3);
            plot(Vlines_red{e}(l).point1(1),Vlines_red{e}(l).point1(2),'ro')
            plot(Vlines_red{e}(l).point2(1),Vlines_red{e}(l).point2(2),'ro')
        end
    end
    
    for e = 1:numel(blue_Vlines)
        for l = 1:numel(blue_Vlines{e})
            L = [blue_Vlines{e}(l).point1;blue_Vlines{e}(l).point2];
            plot(L(:,1),L(:,2),'b','LineWidth',3);
            plot(blue_Vlines{e}(l).point1(1),blue_Vlines{e}(l).point1(2),'ro')
            plot(blue_Vlines{e}(l).point2(1),blue_Vlines{e}(l).point2(2),'ro')
        end
    end
    hold off
    
    
    
    hold off
    imshow(I,[]);
    hold on
    for e = 1:numel(OUT_SQUARE)
        plot(OUT_SQUARE{e}(:,1),OUT_SQUARE{e}(:,2),'g','LineWidth',3);
        plot(IN_SQUARE{e}(:,1),IN_SQUARE{e}(:,2),'g','LineWidth',3);
        plot(RED_SQUARE{e}(:,1),RED_SQUARE{e}(:,2),'r','LineWidth',3);
        plot(BLUE_STRIP_SQUARE{e}(:,1),BLUE_STRIP_SQUARE{e}(:,2),'b','LineWidth',3);
        plot(LEFT_STRIP_SQUARE{e}(:,1),LEFT_STRIP_SQUARE{e}(:,2),'m','LineWidth',3);
        plot(RIGHT_STRIP_SQUARE{e}(:,1),RIGHT_STRIP_SQUARE{e}(:,2),'m','LineWidth',3);
    end
    
    saveas(gca,[oPath num2str(n) '.jpg']);
    pause(.5);
    close all
end























