function [paraData,paraMod,imgX,imgY] = generatePARAcombinations(para,paraSpace,N,I,xD,sz,fsz,pt,disp)
    sel = randperm(size(paraSpace,1));
    sel = sel(1:N);
    paraData = repmat(para,[numel(sel) 1]);
    paraMod = paraData;
    imgX = zeros([fsz N]);

    basis = [[1 0 0]',[0 1 0]',[-1 0 0]',[0 -1 0]'];
    imgY = zeros(N,11);
    for e = 1:numel(sel)
        fprintf(['starting:' num2str(e) '\n']);

        toNaturalSystem = para;
        toNaturalSystem(1) = toNaturalSystem(1);
        toNaturalSystem(2) = toNaturalSystem(2);
        toNaturalSystem(3) = toNaturalSystem(3);
        [naturalT] = generateDomainTransformation([toNaturalSystem(1:3) 0 0]);
        %[naturalT] = toFromaffine(naturalT);
        [naturalD] = transformDomain(xD,naturalT);
        naturalI = myInterp2Sampler(I,pt,naturalD,sz);
        




        paraData(e,1) = -paraSpace(sel(e),1);
        %paraData(e,2) = paraData(e,2)*paraSpace(sel(e),2)^-1;
        %paraData(e,3) = paraData(e,3)*paraSpace(sel(e),3)^-1;
        paraData(e,2) = paraSpace(sel(e),2);
        paraData(e,3) = paraSpace(sel(e),3);


        
        paraMod(e,1) = 0;%paraMod(e,1);
        paraMod(e,2) = (paraSpace(sel(e),2)/paraMod(e,2)).^-1;
        paraMod(e,3) = (paraSpace(sel(e),3)/paraMod(e,3)).^-1;
        [modT] = generateDomainTransformation([paraMod(e,:) 0 0]);
        [xsD] = transformDomain(xD,modT);
        paraMod(e,1) = paraSpace(sel(e),1);
        tmpI = myInterp2Sampler(I,pt,xsD,sz);
        tmpI = imrotate(tmpI,180/pi*paraMod(e,1),'bilinear','crop');
        for r = 1:4
            tmpI(1:15,:) = [];
            tmpI = imrotate(tmpI,90);
        end

        % get the size before zscore normalization
        szr = size(tmpI);
        % zscore normalize the sampled image
        tmpI = zscore(tmpI(:));
        % reshape the transformed image
        tmpI = reshape(tmpI,szr);
       

        %{
        % TRY DISK
        [t1 t2] = ndgrid(linspace(-30,30,61),linspace(-pi,pi,100));
        tD = [t1(:).*cos(t2(:)) t1(:).*sin(t2(:))];
        tryPara = [0 1 1 0 0];
        [tryT] = generateDomainTransformation(tryPara);
        [tD] = transformDomain(tD,tryT);
        [tmpI] = myInterp2Sampler(tmpI,[31 31]',tD,size(t1));
        %}


        imgX(:,:,e) = tmpI;
        [Tx] = generateDomainTransformation([paraData(e,:) 0 0]);
        Ys = Tx*basis;
       

        paraData(e,1)

        fprintf(['ending:' num2str(e) '\n']);
        imgY(e,:) = [Ys(:)' paraData(e,1) paraData(e,2:3)];
        if  disp
            
            imshow(imgX(:,:,e),[]);
            hold on
            plot(Ys(1,:)+(fsz(1)-1)/2+1,Ys(2,:)+(fsz(1)-1)/2+1,'r.')
            plot((fsz(1)-1)/2+1,(fsz(2)-1)/2+1,'r*')

            TH = linspace(-pi,pi,50);
            dispX = [cos(TH)' sin(TH)' ones(50,1)];
            dispX = (Tx*dispX')';
            plot(dispX(:,1)+(fsz(1)-1)/2+1,dispX(:,2)+(fsz(1)-1)/2+1,'g')
            title(num2str(paraData(e,:)))
            drawnow
            hold off
            waitforbuttonpress
        end

    end

    

           
end