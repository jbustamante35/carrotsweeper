function [returnVector,returnPointList,dX,networkOutput,patch] = generatePARAfromNetworks_ver3(I,P,xD,sz,trainedNetworks,Z0,maxIter,disp)

    if numel(P) == 1
        [P(1,2),P(1,1)] = ind2sub(size(I),P);
    end
    displaceMent = [0 0];

    rot = 0;
    

    % produce best guess without grader
    for iter = 1:maxIter

       
        paraSample(1) = mean(rot);
        paraSample(2:3) = 1;
        paraSample(4:5) = 0;
        [Tx] = generateDomainTransformation(paraSample);
        [xDs] = transformDomain(xD,Tx);
        xDs = double(xDs);
        pToSample = P(end,:)';
        [patch] = myInterp2Sampler(I,pToSample,xDs,sz);
        sz = size(patch);
        patch = zscore(patch(:),1,1);
        patch = reshape(patch,sz);


        % apply network to patch
        networkOutput = zeros(1,5);
        for loop = 1:numel(trainedNetworks)
            networkOutput(loop) = trainedNetworks{loop}.predict(patch)*Z0(1,loop) + Z0(2,loop);
        end
        % obtain direct para
        [directPara] = obtainParaFromNetworkOutput(networkOutput);
        % obtain point list
        [indirectPointList] = obtainPointListFromNetworkOutput(networkOutput);
        % obtain a suggested para from point list
        [indirectPara,dX] = obtainParaFromPointList(indirectPointList);
        % obtain directPointList
        [directPointList] = obtainPointListFromPara(directPara);
        % stack the results
        returnVector = [directPara;indirectPara];
        % stack the pointlists
        returnPointList = cat(3,directPointList,indirectPointList);
        


        if disp
            close all
            CL = {'r' 'g' 'b'};
            TYPE = {'.' '*' 'o'};
            imshow(patch,[]);
            hold on
            loc = [31 31];
            for g = 1:size(returnVector,1)
                paraDisplay = returnVector(g,:);
                pointListDisplay = returnPointList(:,:,g);

                displayStoma('',paraDisplay,loc,pointListDisplay',CL{g},TYPE{g},[1 10]);
            end
            title([num2str(iter) '---' num2str(paraDisplay(1)) '---d:' num2str(displaceMent)])
         
           
            %waitforbuttonpress
            quiver(31,31,100*displaceMent(1),100*displaceMent(2),'y')
            drawnow
            hold off
        end
       
       
        rot(iter+1) = rot(iter) + mean(returnVector(:,1),1);
        P(iter+1,:) = P(iter,:) + [1 1].*displaceMent;

    end

    


end