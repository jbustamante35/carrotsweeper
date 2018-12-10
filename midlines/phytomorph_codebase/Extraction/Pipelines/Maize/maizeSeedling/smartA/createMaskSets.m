function [M,LABEL] = createMaskSets(I,GMM1,GMM2)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% find background
    % look at the size of the crop box
    sz = size(I);
    % reshape the data
    tmp1 = reshape(I,[prod(sz(1:2)) sz(3)]);
    [K] = GMM1.cluster(tmp1);
    [K] = reshape(K,sz(1:2));
    BK = K == 1;
    %%%%%%%% find background
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% find foreground
    FG = find(K == 2);
    K2 = GMM2.cluster(tmp1(FG,:));
    M2 = zeros(sz(1:2));
    M2(FG) = K2;
    
    
    %%%%%%%% find foreground
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    M = zeros([sz(1:2) 2]);
    M(:,:,1) = K==1;
    M(:,:,2) = M2==1;
    M(:,:,3) = K==2;
    LABEL = M2;
    %{
    K2 = K2 + 1;

    %%%%%%%% find foreground
    LABEL = zeros(sz(1:2));
    LABEL(find(K==2)) = K2;
    LABEL(find(K==1)) = 1;
    LABEL = label2rgb(LABEL);
               %} 

end