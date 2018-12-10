function [O dis worked] = measureMovie(sm,mm,pcaFlag)
    worked = 1;
    dis = 0;
    try
        firstFrames = 30;
        sz = size(sm);
        %mm = imerode(mm,strel('disk',23,0));
        for e = 1:size(sm,4)
            tmpM = getTopMask(sm(:,:,:,e));
            h(:,:,e) = rgb2hsv_fast(sm(:,:,:,e),'','H');
            sm(:,:,:,e) = bsxfun(@times,sm(:,:,:,e),mm);
            redRatio(e) = sum(tmpM(:)) / sum(mm(:));
            ss(:,:,:,e) = imresize(sm(:,:,:,e),.25);
        end
        szr = size(ss);
        mini = reshape(ss,[prod(szr(1:3)) szr(4)]);



        W = 20;
        x = 1:(W+1)';
        for str = 1:(size(mini,2)-W-1)
            %{
            [mS mC mU mE] = PCA_FIT_FULL_T(mini(:,str:(str+W)),3);
            [mC] = PCA_REPROJ_T(mini(:,str+W+1),mE,mU);
            mC = PCA_BKPROJ_T(mC,mE,mU);
            %}
            %{
            Y = mini(:,str:(str+W))';
            for d = 1:size(Y,2)
                p = polyfit(x',Y(:,d),1);
                val(d) = polyval(p,(x(end)+1));
            end
            %}
            Y = mini(:,str:(str+W))';
            Yn = mini(:,(str+W+1))';
            M = mean(Y(end-4,:),1) - mean(Y(1:5,:),1);
            val = mean(Y(1:5,:),1) + (str+W+1)*M;
            deltaW = (val - Yn);
            deltaP(str) = sum(deltaW.*deltaW,2).^.5;
            %{
            f = fit(x',mini(:,str:(str+W))','polyl');
            f = spap2(1,1,x,mini(:,str:(str+W)));
            val = fnval(f,22);
            b = robustfit(x,mini(:,str:(str+W)));
            p = polyfit(x,mini(:,str:(str+W),1);
            deltaW = (mini(:,str+W+1)-mC);
            deltaP(str) = sum(deltaW.*deltaW,1).^.5;
            %}
        end




        [mS mC mU mE] = PCA_FIT_FULL_T(mini(:,1:firstFrames),3);
        [mC] = PCA_REPROJ_T(mini,mE,mU);
        mC = PCA_BKPROJ_T(mC,mE,mU);
        deltaW = (mini-mC);
        delta = sum(deltaW.*deltaW,1).^.5;
        O = delta;

        deltaP = [delta(1:(numel(delta)-numel(deltaP))) deltaP];

        mU = mU / norm(mU);
        tan = mU'*deltaW;


        O(2,:) = delta - tan;
        O(3,:) = imfilter(redRatio,fspecial('average',[1 11]),'replicate');


        dm = imerode(mm,strel('disk',23,0));
        h = reshape(h,[prod(sz(1:2)) sz(4)]);
        h = h(find(dm==1),:);
        sig = h(:,1);
        u = mean(sig(:)+.5);
        oh = h + (u-.5);
        oh(oh(:) < 0) = 1 +oh(oh(:)<0);
        O(4,:) = mean(oh);
        O(4,:) = abs(subBL(O(4,:)',30)');
        O(5,:) = deltaP;

        %mC = bsxfun(@minus,mC,mU);
        %{
        h = [];
        for e = 1:size(sm,4)
            sm(:,:,:,e) = bsxfun(@times,sm(:,:,:,e),mm);
            h(:,:,e) = rgb2hsv_fast(sm(:,:,:,e),'','H');
        end


        dm = imerode(mm,strel('disk',23,0));
        h = reshape(h,[prod(sz(1:2)) sz(4)]);
        %{
        %m = reshape(m,[prod(sz(1:3)) sz(4)]);
        %mmvec = mm(:);
        h = h(find(dm),:);


        u = mean(h(:)+.5);

        O(1,:) = mean(h,1);
        O(2,:) = mean(h-u,1);
        %{
        plot(O');
        drawnow
        %}
        %}
        h = h(find(dm),:);
        sig = h(:,1:20);
        u = mean(sig(:)+.5);
        oh = h - (u-.5);
        oh(oh(:) < 0) = 1 +oh(oh(:)<0);
        sig = abs(diff(oh,1,2));


        O(1,:) = mean(oh,1);
        O(2,:) = std(oh,1,1);
        O(3,:) = mean(sig,1);
        %{
        sm = reshape(sm,[prod(sz(1:3)) sz(4)]);
        sm = imfilter(sm,fspecial('average',[11 1]),'replicate');
        if pcaFlag
            [S sm U E] = PCA_FIT_FULL_T(sm,5);
        end
        t = gradient(sm);
        n = gradient(t);
        T = bsxfun(@times,t,sum(t.*t,1).^-.5);
        nor = n - bsxfun(@times,sum(n.*T,1),T);
        N = bsxfun(@times,nor,sum(nor.*nor,1).^-.5);
        K = sum(n.*N,1);
        A = sum(n.*T,1);
        L = sum(t.*T,1);
        O = [K;A;L];
        %}

        for e = 1:size(dis,4)
            %imshow(abs(imresize(dis(:,:,:,e),4))*100,[ ]);
            imshow(abs(imresize(sm(:,:,:,e),4)),[ ]);
            title(num2str(e))
            drawnow

        end

        for e = 1:size(sm,4)
            imshow(sm(:,:,:,e),[]);
            title(num2str(e))
            drawnow

        end
        %}
    catch 
        worked = 0;
    end
end