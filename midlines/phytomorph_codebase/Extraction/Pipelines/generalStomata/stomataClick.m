function [xData,yData,selP,pVEC,IDXMASTER,PARAMASTER] = stomataClick(result,Grade,iFileList,chop)

    xData = [];
    yData = [];
    
    for e = 1:numel(result)
        IDXMASTER{e} = [];
        PARAMASTER{e} = [];


        % read the image
        testRAW = imread(iFileList{e});
        test = testRAW;

        % snip the borders
        for r = 1:4
            test(1:chop,:) = [];
            test = imrotate(test,90);
        end



        TESTER = squeeze(result{e}(:,end-3:end-1));
        TESTER = reshape(TESTER,[size(test) size(TESTER,2)]);

        G = Grade{e};
        G = [prod(G(:,1:2),2) prod(G(:,3:4),2) prod(G(:,5:6),2)];

        CL = {'r' 'g' 'b'};
        figure;
        imshow(test,'Border','tight','DisplayRange',[]);
        hold on
        totalPoints = 0;
        hx1T = [];
        hy1T = [];
        for g = 1:size(G,2)
            
            tG = reshape(G(:,g),size(test));
            hope = -log(tG);
            msk = imerode(hope,strel('disk',11,0))==hope;
            msk = msk==1 & bwareaopen(msk,3)==0;
            [hx1,hy1] = find(msk);
           
            hy1T = [hy1T;hy1];
            hx1T = [hx1T;hx1];


            plot(hy1,hx1,'.');

            for p = 1:numel(hx1)
                para = squeeze(TESTER(hx1(p),hy1(p),:));
                para = [para;0;0];
                [T] = generateDomainTransformation(para);
                T(:,3) = T(:,3) + [hy1(p);hx1(p)];
                TH = linspace(-pi,pi,100);
                dX10 = [cos(TH)' sin(TH)'];
                [dX1] = transformDomain(dX10,T);
                plot(dX1(:,1),dX1(:,2),[CL{g} '--']);

                IDX = sub2ind(size(test),hx1(p),hy1(p));




                tmp = result{e}(IDX,1:end-4);
                xData = [xData;tmp];

            end
            totalPoints = totalPoints + numel(hy1);
        end








        selP{e} = [];
        [selP{e}(:,1),selP{e}(:,2),V] = impixel();
        pVEC{e} = zeros(totalPoints,1);
        rThresh = 3;
        % for each selected point - find the nearest and those less than r
        for p = 1:size(selP{e},1)

            delta = bsxfun(@minus,[hy1T hx1T],selP{e}(p,:));
            delta = sum(delta.*delta,2).^.5;
            midx = find(delta < rThresh);
            %[~,midx] = min(delta);

            % set value to true for those within rThresh
            pVEC{e}(midx) = 1;

            for m = 1:numel(midx)
                para = squeeze(TESTER(hx1T(midx(m)),hy1T(midx(m)),:));
                para = [para;0;0];
                [T] = generateDomainTransformation(para);
                T(:,3) = T(:,3) + [hy1T(midx(m));hx1T(midx(m))];
                TH = linspace(-pi,pi,100);
                dX10 = [cos(TH)' sin(TH)'];
                [dX1] = transformDomain(dX10,T);
                plot(dX1(:,1),dX1(:,2),'g');


                

                IDXMASTER{e} = [IDXMASTER{e};sub2ind(size(testRAW),hx1T(midx(m))+chop,hy1T(midx(m))+chop)];
                PARAMASTER{e} = [PARAMASTER{e};para];
            end





        end
        drawnow
        waitforbuttonpress
        yData = [yData;pVEC{e}];
       
        close all
    end
end