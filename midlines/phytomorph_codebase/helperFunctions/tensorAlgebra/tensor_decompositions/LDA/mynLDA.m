function [lambda] = mynLDA(D,G,curL,threshL)
    C = 0;
    UNI = unique(G);
    for i = 1:2
        tD = D(G==UNI(i),:);
        U(i,:) = mean(tD,1);
        tD = bsxfun(@minus,tD,U(i,:));
        %tD = tD - repmat(U(i,:),[size(tD,1) 1]);
        %C = C + tD'*tD;
        C = C + mtimesx(tD,'T',tD);
    end
    delta = diff(U,1,1)';    
    lambda = C\delta;
    lambda = lambda/norm(lambda);
    if curL == threshL
        return
    else
        [D,BV] = myGS(D,lambda);        
        lambdaN = mynLDA(D,G,curL+1,threshL);
        lambda = [lambda BV*lambdaN];
    end
end
%{
    load fisheriris
    testD = meas(strcmp(species,'setosa') | strcmp(species,'versicolor'),:);
    testG = zeros(size(testD,1),1);
    idx1 = find(strcmp(species,'setosa'));
    idx2 = find(strcmp(species,'versicolor'));
    testG(idx1) = 1;
    testG(idx2) = 2;
    lambda = myLDA(testD,testG);

%}


%http://digital.library.adelaide.edu.au/dspace/bitstream/2440/15227/1/138.pdf
