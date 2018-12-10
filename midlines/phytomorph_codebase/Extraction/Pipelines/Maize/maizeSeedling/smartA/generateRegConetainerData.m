function [modelTable] = generateRegConetainerData(masterTable,modelTable)

    %% find masked conetainer data
    midx = find(strcmp(masterTable.type,'conetainer_whole_masked'));
    basePath = '/mnt/tetra/nate/seedlingImageParts/conetainers/registeredModels/';
    ZS = [];
    cnt = 1;
    topColor = [];
    for e = 1:numel(midx)
        fileName = masterTable.imageLocation{midx(e)};
        [pth,nm,ext] = fileparts(fileName);
        I = imread(fileName);

        wholeFile = [pth filesep nm ext];
        wholeFile = strrep(wholeFile,'_masked','');
        wI = imread(wholeFile);
        

        msk = all(I ~= 0,3);
        msk = imfill(msk,'holes');
        BW = edge(msk);
        BW = imdilate(BW,strel('disk',3,0));

        TH = linspace(-15,15,100);
        [H,T,R] = hough(BW','Theta',TH);
        P  = houghpeaks(H,1);
        lines = houghlines(BW',T,R,P,'FillGap',1000,'MinLength',.5*size(I,2));

        
        if ~isempty(lines)
            k = 1;
            xy = [lines(k).point1; lines(k).point2];
            slope = diff(xy,1,1);
            angle = atan2(slope(1),slope(2));
            oldANGLE = angle*180/pi;
            I = imrotate(I,angle*180/pi);
            wI = imrotate(wI,angle*180/pi);
            BW = imrotate(BW,angle*180/pi);
            msk = imrotate(msk,angle*180/pi);
            TH = linspace(-15,15,1000);
            [H,T,R] = hough(BW','Theta',TH);
            P  = houghpeaks(H,1);
            lines = houghlines(BW',T,R,P,'FillGap',1000,'MinLength',.5*size(I,2));
            if ~isempty(lines)
                
                
                xy = [lines(k).point1; lines(k).point2];
                slope = diff(xy,1,1);
                angle = atan2(slope(1),slope(2));
                imshow(I,[]);  hold on
                m = slope(1)/slope(2);
                
                midpoint = mean(xy,1);
                
                b = midpoint(1) - m*midpoint(2);
                toPlot = [[1 b];[size(I,2) b + size(I,2)*m]];
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                coneMask = msk;
                coneMask(round(b)-5:end,:) = 0;
                coneMask = bwlarge(coneMask);
                WID = sum(coneMask,2);
                WID1 = median(WID(WID~=0));
                WID2 = mean(WID(WID~=0));
                WID = .5*(WID1 + WID2);
                R = regionprops(logical(coneMask),'Centroid');
                coneBOX = [R(1).Centroid(1)-WID/2 1 WID b];
                rectangle('Position',coneBOX,'EdgeColor','r')
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                G = rgb2gray(I);
                bE = edge(G,'Canny');
                LB = (R(1).Centroid(1)-WID/2 + WID/2) - 1.5*WID/2;
                RB = (R(1).Centroid(1)-WID/2 + WID/2) + 1.5*WID/2;
                subG = G(:,round(LB:RB));
                [dG1 dG2] = gradient(double(subG));
                
                dG1(1:b,:) = 0;
                dG2(1:b,:) = 0;
                %{
                plot(sum(dG1,2))
                hold all
                plot(sum(dG2,2))
                figure;
                %}
                
                if e < 5
                    
                    [c r V] = impixel(G);
                    subG = imresize(subG,[size(subG,1) 301]);
                    tmpF = subG(r(1):r(2),:);
                    tmpF = imresize(tmpF,[100 301]);
                    sigS(:,:,e) = double(tmpF);
                    
                    fun = sum(dG1,2).^-1 .* sum(dG2,2);
                    [~,loc] = max(fun);
                    LL = [];
                    RR = [];
                    
                    
                    
                else
                    
                    
                    fSIG = mean(sigS,3);
                    fSIG = fSIG(:);
                    fSIG = fSIG / norm(fSIG);
                    subG = imresize(subG,[size(subG,1) 301]);
                    BL = double(im2col(subG,[100 301],'sliding'));
                    n = sum(BL.*BL,1).^-.5;
                    BL = bsxfun(@times,BL,n);
                    fSIG = fSIG'*BL;
                    [~,loc] = max(fSIG);
                    yO = mean(toPlot(:,2));
                    smsk = msk(yO-2:yO+2,:);
                    msig = sum(smsk,1);
                    msig = msig > 3;
                    msig = [ones(size(msig));msig;ones(size(msig))];
                    msig = imfill(msig,'holes');
                    msig = msig(2,:);
                    msig = find(msig);
                    
                    
                    
                    
                    otherSIG = [[1 loc];[size(I,2) loc]];
                    
                    
                    outV = imref2d([1000,1000]);
                    
                    
                    
                    imshow(G,[]);
                    hold on;
                    plot(toPlot(:,1),toPlot(:,2),'r');
                    plot(otherSIG(:,1),otherSIG(:,2),'r');
                    [c r V] = impixel();
                    
                    
                    movingPoints = fliplr([r c]);
                    fixedPoints = fliplr([[1 1];[1 1000];[1000 1];[1000 1000]]);
                    tform = fitgeotrans(movingPoints,fixedPoints,'projective');
                    
                    tmpP = G(b:loc,:);
                    tmpP = imresize(tmpP,[100,500]);
                    toPredict(:,:,cnt) = tmpP;
                    
                    TF{cnt} = tform;
                    topColor(:,:,:,cnt) = imwarp(wI,tform,'OutputView',outV);
                    cnt = cnt + 1;
                    
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                otherSIG = [[1 loc];[size(I,2) loc]];
             
                
                
                
                %{
                eC = edge(coneMask);
                eC = imdilate(eC,strel('square',3));
                TH = linspace(-15,15,100);
                [H,T,R] = hough(eC,'Theta',TH);
                P = houghpeaks(H,10,'NHoodSize',[101 21]);
                lines = houghlines(eC,T,R,P,'FillGap',100);
                for k = 1:numel(lines)
                    l1 = [lines(k).point1; lines(k).point2];
                    plot(l1(:,1),l1(:,2),'c');
                end
                %}
                
              
                plot(toPlot(:,1),toPlot(:,2),'r');
                if ~isempty(LL)
                    plot(LL(:,1),LL(:,2),'r');
                    plot(RR(:,1),RR(:,2),'r');
                end
                plot(otherSIG(:,1),otherSIG(:,2),'r');
                title([num2str(angle*180/pi) '--' num2str(oldANGLE)]);
                drawnow
                
                
                hold off
                pause(.21);
                if angle == 0
                    stop = 1;
                end
                pt = size(modelTable,1);
                modelTable{pt+1,'type'} = {'registeredConetainer'};
                modelTable{pt+1,'imageLocation'} = {[basePath nm '.tif']};
                for k = 1:4
                    newKey =  ['original_boundingBox' num2str(k)];
                    oldKey = ['boundingBox' num2str(k)];
                    modelTable{pt+1,newKey} = masterTable{midx(k),oldKey};
                end
                modelTable{pt+1,'modelAngle'} = angle;
                modelTable{pt+1,'modelWidth'} = size(I,2);
                modelTable{pt+1,'modelHeight'} = size(I,1);
                waitforbuttonpress
                
            end
        end
        
        
    end
end