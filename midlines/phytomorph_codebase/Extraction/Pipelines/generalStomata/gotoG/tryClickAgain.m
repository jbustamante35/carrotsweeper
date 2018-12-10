function [] = tryClickAgain(oI,fI)

    h1 = figure;
    h2 = figure;
    toclick = true;
    
    X = [];
    Y = [];
    
    for e = 1:size(oI,3)


        pleaseWORK = fI(:,:,:,e);
        pleaseWORK = reshape(pleaseWORK,[size(pleaseWORK,1)*size(pleaseWORK,2) size(pleaseWORK,3)]);
        
        if e > 2
            try
                %[labels, posterior, cost]
                [TYPE,PROB] = predict(WHY,pleaseWORK);
            catch ME
                ME = 1;
            end
            
            PROB = sum(PROB(:,1:cluster1),2);

            %{
            PROB = mvnpdf(pleaseWORK(:,toUse),mean(SC(:,toUse)),cov(SC(:,toUse)));
            %{
            PROB = bsxfun(@minus,pleaseWORK,mean(SC,1));
            PROB = bsxfun(@mtimes,PROB,std(SC,1,1).^-1);
            PROB = (sum(PROB.*PROB,2)).^.5;
            %}
            %}
            TYPE = reshape(TYPE,[size(oI,1) size(oI,2)]);
            PROB = reshape(PROB,[size(oI,1) size(oI,2)]);
            
            figure(h1);
            imshow([bindVec(PROB) bindVec(oI(:,:,e))],[]);
            %{
            PROB = log(.01*PROB);
            PROB(isinf(PROB)) = 0;
            PROB = imfilter(PROB,fspecial('gaussian',[31 31],5),'replicate');
            PROB = bindVec(PROB);


            if ~isempty(PROB)
                figure(h1);
                imshow(PROB,[]);
                if ~toclick
                    figure(h2)
                    out = flattenMaskOverlay(bindVec(oI(:,:,e)),reshape(TYPE,[size(oI,1) size(oI,2)])<= cluster1);
                    imshow(out,[]);
                    waitforbuttonpress
                end

            end
            %}
        end


        if toclick
            figure(h2)
            [cp{e}(:,2),cp{e}(:,1),~] = impixel(oI(:,:,e),[]);
            
            
            
            IDX = sub2ind([size(oI,1) size(oI,2)],cp{e}(:,1),cp{e}(:,2));
            msk = zeros(size(oI,1),size(oI,2));
            msk(IDX) = 1;
            msk = imdilate(msk,strel('disk',5,0));
            IDX = find(msk==1);
            out = flattenMaskOverlay(bindVec(oI(:,:,e)),logical(msk));
            imshow(out,[]);
            drawnow
            waitforbuttonpress
            
            %{
            FUN = [];
            for p = 1:size(cp{e},1)
                FUN= cat(3,FUN,squeeze(retVec1(cp{e}(p,2),cp{e}(p,1),:,:,2,e)));
            end
            %}
            
            Y = [Y;msk(:)];
            X = [X;pleaseWORK];


            if e >=2
               

                F1 = Y==1;
                
                %{
                F1 = Y==1;
                AIC = [];
                BIC = [];
                GMModel_sub = {};
                for cluster1 = 1:4
                    GMModel_sub{cluster1} = fitgmdist(X(F1,:),cluster1,'Options',options,'Replicates',3);
                    AIC(cluster1) = GMModel_sub{cluster1}.AIC;
                    BIC(cluster1) = GMModel_sub{cluster1}.BIC;
                end
                [minAIC,cluster1] = min(AIC);
                %}
                
                cluster1 = 3;
                clear GMModel_sub;
                
                options = statset('Display','iter','MaxIter',300);
                GMModel_sub = fitgmdist(X(F1,:),cluster1,'Options',options,'Replicates',3);
                [subIDX1] = GMModel_sub.cluster(X(F1,:));

                
                
                
                
                
                cluster2 = 3;
                F0 = find(Y==0);
                sam = 10;
                GMModel_sub = fitgmdist(X(F0(1:sam:end),:),cluster2,'Options',options,'Replicates',3);
                [subIDX2] = GMModel_sub.cluster(X(Y==0,:));

                YT = Y;
                YT(F1) = subIDX1;
                YT(F0) = subIDX2 + cluster1;
                classNames = {'s','n'};
                WHY = fitcnb(X,YT);


            end
        end
    end
end