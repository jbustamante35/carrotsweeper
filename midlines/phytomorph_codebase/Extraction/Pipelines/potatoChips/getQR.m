function [qr,msg] = getQR(I,B)
    G = rgb2gray(I);
    E = edge(G,'canny');
    E(1:100,1:100) = 0;
    
    Ef = bwareaopen(E,300);
    E = xor(Ef,E);
    E = (E - B) == 1;
    MASK = zeros(size(E));
    MASK(1:1500,1:1500) = 1;
    E = E.*MASK;
    %{
    R = regionprops(logical(E),'Solidity','PixelIdxList');
    RAT = ([R.Solidity] > .8) & ([R.Solidity] < 1.2);
    R = R(RAT);
    E = zeros(size(E));
    for e = 1:numel(R)
        E(R(e).PixelIdxList) = 1;
    end
    %}
    
    
    
    E = imclose(E,strel('disk',3));
    Eh = imfill(E,'holes');
    Eh = imopen(Eh,strel('disk',11));
    Eh = bwareaopen(Eh,50^2);
    R = regionprops(Eh,'Perimeter','Area','ConvexArea','PixelIdxList','Centroid');
    N = [];
    for e = 1:numel(R)
        N(e) = norm(R(e).Centroid);
    end
    rm = find(N > 1500);
    R(rm) = [];
    RAT = (abs([R.ConvexArea].*[R.Area].^-1)).*abs(([R.Area] - 70^2));
    [SRAT iRAT] = sort(RAT);
    
    for e = 1:3
        BOXES(e,:) = R(iRAT(e)).Centroid;
    end
    imshow(G,[]);
    hold on
    plot(BOXES(:,1),BOXES(:,2),'r*')
    [J,idx]  = min(sum(BOXES.*BOXES,2));
    UL = BOXES(idx,:);
    UL = UL - 70;
    CB = [UL 300 300];
    qr = imcrop(G,CB); 
    msg = decode_qr(qr);
end