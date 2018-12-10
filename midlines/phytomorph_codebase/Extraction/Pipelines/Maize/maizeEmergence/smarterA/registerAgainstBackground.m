function [cameraShift] = registerAgainstBackground(FileList,rowCount,deltaThresh,tform,uCL,covCL)
        % read the first image
        %I1 = imread(FileList{2});
        I1 = getRectifiedImage(FileList{2},tform)*255;
        Iend = getRectifiedImage(FileList{end},tform)*255;


        %myReg(I1/255,Iend/255,uCL,covCL)



        HSV1 = rgb2hsv(I1);    
        HSVend = rgb2hsv(Iend);    



        % convert to gray scale
        G1 = 255*rgb2gray(I1/255);
        Gend = 255*rgb2gray(Iend/255);
        % get the background
        S1 = G1(1:rowCount,:);
        Send = Gend(1:rowCount,:);
        

        % 
        M1 = S1 > 50;%255*graythresh(S1);
        mM1 = bwareaopen(M1,200);
        sm1 = sum(mM1,1);
        sm2 = sum(mM1,2);
        th1 = sm1 > max(sm1)*.950;
        th2 = sm2 > max(sm2)*.95;
        th1 = bwlarge(th1);
        th2 = bwlarge(th2);
        th = double(th2)*double(th1);
        %dist = bwdist(~M1);
        %th = dist > 70;
        R = regionprops(th,'BoundingBox');

        [rec] = findMaxRec(mM1);
        %{
        % make the mask
        M1 = S1 > 20;%255*graythresh(S1);
        Mend = Send > 20;%255*graythresh(Send);
        % make the final mask
        Mtotal = sum(M1,1) > 20 & sum(Mend,1) > 20;
        %Mtotal = all(M1,1).*all(Mend,1);
        Mtotal = imclose(Mtotal,strel(ones(1,200)));
        Mtotal = bwlarge(Mtotal);
        % get the section to register
        F1 = S1(:,find(Mtotal));
        Fend = Send(:,find(Mtotal));
        % configure the registration
       
        %}
        [optimizer, metric] = imregconfig('monomodal');
        optimizer.MaximumIterations = 800;
        optimizer.MinimumStepLength = 5e-4;
        F1 = imcrop(S1,rec);
        Fend = imcrop(Send,rec);

        szH = size(I1)/2;

        CP = [(rec(1)+(rec(1)+rec(3)))/2 (rec(2)+(rec(2)+rec(4)))/2];
        CPT = eye(3);
        CPT(3,1:2) = -(CP) + [szH(2) szH(1)];
        DIS = affine2d(CPT);

        CP = [(rec(1)+(rec(1)+rec(3)))/2 (rec(2)+(rec(2)+rec(4)))/2];
        CPT = eye(3);
        CPT(3,1:2) = (CP) - [szH(2) szH(1)];
        iDIS = affine2d(CPT);
        % register the first and last
        tformTest = imregtform(Fend, F1, 'affine', optimizer, metric,'DisplayOptimization',true);
        %{
        FT = (DIS.T*tformTest.T);
        tformTest2 = affine2d(FT);
        sz = size(I1);
        for k = 1:3
            ret2(:,:,k) = imwarp(Iend(:,:,k),tformTest2,'OutputView',imref2d(sz(1:2)));
        end

        sz = size(I1);
        for k = 1:3
            ret1(:,:,k) = imwarp(I1(:,:,k),DIS,'OutputView',imref2d(sz(1:2)));
        end

        FF = imwarp(Fend,tformTest,'OutputView',imref2d(size(F1)));

        %}
        % measure the total displacement
        delta = norm(tformTest.T(end,1:2));
        % if total displacement is greater than threshold


        if delta > deltaThresh
            % register each image via parfor
            parfor e = 2:numel(FileList)
                tmpI = getRectifiedImage(FileList{e},tform);
                % convert to gray and register to the second image
                tmpS = rgb2gray(tmpI(1:rowCount,:,:));
                tmpS = imcrop(tmpS,rec);

                [output] = dftregistration(fft2(F1),fft2(tmpS*255),4);
                T = [eye(3)];
                T(3,1) = output(4);
                T(3,2) = output(3);
                cameraShift{e} = affine2d(T);
                % store the camera shift
                %cameraShift{e} = imregcorr(tmpS*255, F1, 'rigid');
                %{
                sz = size(Iend);
                retI = [];
                for k = 1:3
                    retI(:,:,k) = imwarp(Iend(:,:,k),cameraShift{e},'OutputView',imref2d(sz(1:2)));
                end
                %}

            end
        else
            cameraShift = [];
        end
end