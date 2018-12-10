function [] = scanalyzerMain(FileList,tform,tform2,oPath,rPath)
    try
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % load and warp the images
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fidx = find(contains(FileList,[filesep 'vis' filesep]));
        FileList = FileList(fidx);
        fFileList = strrep(FileList,[filesep 'vis' filesep],[filesep 'fluo' filesep]);
        nFileList = strrep(FileList,[filesep 'vis' filesep],[filesep 'nir' filesep]);






        %{
        for e = 1:numel(FileList)
            FileList{e}
            fFileList{e}
            nFileList{e}
        end
        %}


        mkdir(oPath);

        for e = 1:numel(FileList)
            try
                exportFileList = {};
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % load and warp the images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:extracting file name(s) and creating folders. \n']);tic;
                % get file information
                [pth,nm,ext] = fileparts(FileList{e});
                sidx = strfind(pth,filesep);
                %pth = pth(sidx(end)+1:end);
                if isdeployed()
                    pth = pth(1:sidx(1)-1);
                else
                    pth = pth(sidx(end-1)+1:sidx(end)-1);
                end
                TMPoPath = [oPath pth filesep];
                mkdir(TMPoPath);
                RGBpath = [TMPoPath 'RGB' filesep];
                FLUOpath = [TMPoPath 'FLUO' filesep];
                IRpath = [TMPoPath 'IR' filesep];
                mkdir(RGBpath);
                mkdir(FLUOpath);
                mkdir(IRpath);
                nm = [pth '__' nm];
                fprintf(['end:extracting file name(s) and creating folders.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % load and warp the images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % load and warp the images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:loading and registering the images.\n']);tic;
                I = double(imread(FileList{e}))/255;
                I2 = double(imread(fFileList{e}))/255;
                I3 = double(imread(nFileList{e}))/255;
                R = imref2d(size(I));
                I2 = imwarp(I2,tform,'OutputView',R);
                I3 = imwarp(I3,tform2,'OutputView',R);
                fprintf(['end:loading and registering the images.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % load and warp the images
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the mask and contour from the RGB image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:extracting mask.\n']);tic;
                HSV = rgb2hsv(I);
                G = HSV(:,:,2);
                M = G > graythresh(G);
                M = imclearborder(M);
                M = imfill(M,'holes');
                M = bwlarge(M);
                M = imclose(M,strel('disk',41,0));
                M = imopen(M,strel('disk',51,0));
                M = bwlarge(M);
                R = regionprops('table',M,'Area','BoundingBox','Perimeter');
                fprintf(['end:extracting mask.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % make the mask and contour from the RGB image
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % extract and smooth the contour
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:extracting and smoothing contours.\n']);tic;
                dB = bwboundaries(M);
                mI = bsxfun(@times,I,double(M));
                ker_dB{e} = dB{1};
                sm_dB{e} = imfilter(ker_dB{e},fspecial('average',[101 1]),'circular');
                fprintf(['end:extracting and smoothing contours.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % extract and smooth the contour
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the geometry of the kernel
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:extracting kernel geometry.\n']);tic;
                [mm1,mm2,mm3,idx1,idx2] = findKernelTip_chop(sm_dB{e},600,60,mI,false);

                midRib = [mm2(2,:);mm3(2,:)];
                vmidRib = midRib;
                NOR = diff(midRib,1,1);
                kernelMajor = norm(NOR);
                NOR = [NOR(2) -NOR(1)];
                NOR = NOR / norm(NOR);
                midRib = arcLength(midRib,'arcLen',1);

                NOR1_p = (1:2000)' * NOR;
                NOR1_n = -NOR1_p;
                fprintf(['end:extracting kernel geometry.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % find the geometry of the kernel
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SAMPLE IMAGES IN SQUARE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:samping images.\n']);tic;
                plot(sm_dB{e}(:,2),sm_dB{e}(:,1),'r');
                hold on;
                K1 = [];
                K2 = [];
                K3 = [];
                G = [];
                WIDTH = [];
                NORT = {};
                for m = 1:size(midRib,1)

                    NOR1_p_tmp = bsxfun(@plus,flipdim(NOR1_p,1),midRib(m,:));
                    NOR1_n_tmp = bsxfun(@plus,NOR1_n,midRib(m,:));


                    [D,IDX] = pdist2(sm_dB{e},NOR1_p_tmp,'euclidean','Smallest',1);
                    [~,s] = min(D);
                    NOR1_p_tmp(1:s-1,:) = [];

                    [D,IDX] = pdist2(sm_dB{e},NOR1_n_tmp,'euclidean','Smallest',1);
                    [~,s] = min(D);
                    NOR1_n_tmp(s+1:end,:) = [];

                    NORT{m} = [NOR1_p_tmp;midRib(m,:);NOR1_n_tmp];
                    WIDTH(m) = norm(NORT{m}(end,:) - NORT{m}(1,:));
                    NORQ = arcLength(NORT{m},'spec',300);

                    %{
                    plot(NORT{m}(:,2),NORT{m}(:,1),'b')
                    drawnow
                    %}

                    K1(:,m,:) = ba_interp2(I,NORQ(:,2),NORQ(:,1));
                    K2(:,m,:) = ba_interp2(I2,NORQ(:,2),NORQ(:,1));
                    K3(:,m,:) = ba_interp2(I3,NORQ(:,2),NORQ(:,1));
                    G(:,m,1) = NORQ(:,1);
                    G(:,m,2) = NORQ(:,2);
                    %{
                    imshow(K,[]);
                    drawnow
                    %}
                end
                WID = .5*WIDTH;
                k1 = bsxfun(@times,squeeze(sum(K1,1)),WIDTH'.^-1);
                k2 = bsxfun(@times,squeeze(sum(K2,1)),WIDTH'.^-1);
                k3 = bsxfun(@times,squeeze(sum(K3,1)),WIDTH'.^-1);
                k1(WIDTH < 3,:) = 0;
                k2(WIDTH < 3,:) = 0;
                k3(WIDTH < 3,:) = 0;
                exportFileList{end+1} = [RGBpath nm '_SQUARE.tif'];
                imwrite(K1,exportFileList{end});
                exportFileList{end+1} = [FLUOpath nm '_SQUARE.tif'];
                imwrite(K2,exportFileList{end});
                exportFileList{end+1} = [IRpath nm '_SQUARE.tif'];
                imwrite(K3,exportFileList{end});
                exportFileList{end+1} = [TMPoPath nm '_extra.mat'];
                save(exportFileList{end},'G')
                fprintf(['end:samping images.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % SAMPLE IMAGES IN SQUARE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:RGB image display.\n']);tic;
                SKIP = 50;
                image(I);
                %imshow(I,[],'Border','tight');
                hold on
                plot(vmidRib(:,2),vmidRib(:,1),'b');
                plot(sm_dB{e}(:,2),sm_dB{e}(:,1),'r','LineWidth',1);
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'r*')
                %{
                plot(mm1(:,2),mm1(:,1),'k');
                plot(mm2(:,2),mm2(:,1),'b');
                plot(mm3(:,2),mm3(:,1),'c');

                plot(sm_dB{e}(idx1,2),sm_dB{e}(idx1,1),'b*')
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'c*')
                %}
                for level = 1:SKIP:numel(NORT)
                    plot(NORT{level}(:,2),NORT{level}(:,1),'r')
                end
                hold off

                drawnow

                exportFileList{end+1} = [RGBpath nm '.jpg'];
                saveas(gca,exportFileList{end});
                close all
                fprintf(['end:RGB image display.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:FLUO image display.\n']);tic;
                SKIP = 50;
                image(I2);
                %imshow(I,[],'Border','tight');
                hold on
                plot(vmidRib(:,2),vmidRib(:,1),'b');
                plot(sm_dB{e}(:,2),sm_dB{e}(:,1),'r','LineWidth',1);
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'r*')
                %{
                plot(mm1(:,2),mm1(:,1),'k');
                plot(mm2(:,2),mm2(:,1),'b');
                plot(mm3(:,2),mm3(:,1),'c');

                plot(sm_dB{e}(idx1,2),sm_dB{e}(idx1,1),'b*')
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'c*')
                %}
                for level = 1:SKIP:numel(NORT)
                    plot(NORT{level}(:,2),NORT{level}(:,1),'r')
                end
                hold off

                drawnow

                exportFileList{end+1} = [FLUOpath nm '.jpg'];
                saveas(gca,exportFileList{end});
                close all
                fprintf(['end:FLUO image display.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:NIR image display.\n']);tic;
                SKIP = 50;
                image(I3);
                %imshow(I,[],'Border','tight');
                hold on
                plot(vmidRib(:,2),vmidRib(:,1),'b');
                plot(sm_dB{e}(:,2),sm_dB{e}(:,1),'r','LineWidth',1);
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'r*')
                %{
                plot(mm1(:,2),mm1(:,1),'k');
                plot(mm2(:,2),mm2(:,1),'b');
                plot(mm3(:,2),mm3(:,1),'c');

                plot(sm_dB{e}(idx1,2),sm_dB{e}(idx1,1),'b*')
                plot(sm_dB{e}(idx2,2),sm_dB{e}(idx2,1),'c*')
                %}
                for level = 1:SKIP:numel(NORT)
                    plot(NORT{level}(:,2),NORT{level}(:,1),'r')
                end
                hold off
                drawnow
                exportFileList{end+1} = [IRpath nm '.jpg'];
                saveas(gca,exportFileList{end});
                close all
                fprintf(['end:NIR image display.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % DISPLAY THE IMAGE
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % square-interpolate the 
                icon = circshift(sm_dB{e},-(idx2-1),1);
                iContour = interp1(1:size(sm_dB{e},1),icon,linspace(1,size(sm_dB{e},1),2000));

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %% spin-up JSON output format
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:JSON export.\n']);tic;
                phDoc = [];
                phDoc = generatePhenotypeNode(phDoc,[pth filesep nm ext],'RGB_fileName','RGB_fileName');
                phDoc = generatePhenotypeNode(phDoc,[pth filesep nm ext],'FLUO_fileName','FLUO_fileName');
                phDoc = generatePhenotypeNode(phDoc,[pth filesep nm ext],'NIR_fileName','NIR_fileName');
                phDoc = generatePhenotypeNode(phDoc,R{1,'Area'},'kernelArea','kernelArea');
                phDoc = generatePhenotypeNode(phDoc,R{1,'Perimeter'},'kernelPerimeter','kernelPerimeter');
                phDoc = generatePhenotypeNode(phDoc,kernelMajor,'kernelLength','kernelLength');
                phDoc = generatePhenotypeNode(phDoc,WIDTH,'kernelWidthProfile','kernelWidthProfile');
                phDoc = generatePhenotypeNode(phDoc,k1,{'along_kernelAxis' 'colorDims'},'RGBprofile');
                phDoc = generatePhenotypeNode(phDoc,k2,{'along_kernelAxis' 'colorDims'},'FLUOprofile');
                phDoc = generatePhenotypeNode(phDoc,k3,{'along_kernelAxis' 'colorDims'},'NIRprofile');
                phDoc = generatePhenotypeNode(phDoc,iContour,{'around_kernelPerimeter' 'x_yDim'},'interpolatedContour');
                JSON_string = savejson('kernelDoc',phDoc);


                %%%%%%%%%%%%%%%%%%%%%%%%%
                % save JSON string
                %%%%%%%%%%%%%%%%%%%%%%%%%
                exportFileList{end+1} = [TMPoPath filesep nm '_jdoc.json'];
                fileID = fopen(exportFileList{end},'w');
                fprintf(fileID,strrep(JSON_string,'\/','\\/'));
                fclose(fileID);
                %%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['end:JSON export.(' num2str(toc) ':secs)\n'])

                %%%%%%%%%%%%%%%%%%%%%%%%%
                % push to iRODS
                %%%%%%%%%%%%%%%%%%%%%%%%%
                fprintf(['start:IRODS export.\n']);tic;
                pushToiRods(rPath,exportFileList);
                fprintf(['end:IRODS export.(' num2str(toc) ':secs)\n'])
                %%%%%%%%%%%%%%%%%%%%%%%%%
            catch
            end
        end



         %{
         [tipPoint,dB] = getInitialGuessForTip_scanalyzer(sm_dB);
         for e = 1:numel(FileList)
            I = double(imread(FileList{e}))/255;
            imshow(I,[]);
            hold on
            plot(sm_dB{e}(:,2),sm_dB{e}(:,1),'r','LineWidth',2);
            plot(tipPoint(e,2),tipPoint(e,1),'r*')
            title(num2str(e));
            drawnow
            waitforbuttonpress

         end
         %}
    catch ME
        getReport(ME)
    end
end


%{
    func = @(X)scanalyzerMain(X,tform,tform2,'./output/','');
    f = partialFunction(func,'scanalyzer');
    f.publish();
%}

