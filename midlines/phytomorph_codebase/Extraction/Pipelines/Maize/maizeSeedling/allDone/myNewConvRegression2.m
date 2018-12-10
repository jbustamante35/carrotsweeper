function [D] = myNewConvRegression2(K,X,Y,SZ,xSZ,P,h,E,U)

    %K = K / sum(K);
    %K = K / norm(K);
    CONV = mtimesx(K',X);

    for tr = 1:xSZ(4)
       
        TCONV = col2im(CONV(:,:,tr),SZ(1:2),[xSZ(1:2)]);
        %TCONV = 1 + TCONV;
        TCONV(TCONV<0) = 0;
        TCONV = TCONV * sum(sum(TCONV,1),2)^-1;
        
        Q(tr,1) = sum(sum(P(:,:,1).*TCONV,1),2);
        Q(tr,2) = sum(sum(P(:,:,2).*TCONV,1),2);
        
        
        if tr == 1
            %K = PCA_BKPROJ_T(K,E,U);
            %figure(h)
            %imshow(reshape(K,[14 14 3]),[])
            
            %imshow(TCONV,[]);
            %imshow(K,[]);
            %drawnow
        end
        
        
        
    end
    D = mean(sum((Q - Y).^2,2).^.5,1);
end

%{
    P = [];
    SZ = [14 14 3];
    HW = (SZ-1)/2;
    [P(:,:,1),P(:,:,2)] = ndgrid(HW(1)+1:(size(NEWX,1)-HW(1)),HW(2)+1:(size(NEWX,2)-HW(2)));
    
    options = optimset('Display','iter');

    x0 = ones(prod(SZ),1)*prod(SZ)^-1;


    X = [];
    SZ = [14 14 3];
    for tr = 1:size(NEWX,4)
        tmp = NEWX(:,:,:,tr);
        tmp = im2colF(tmp,SZ,[1 1 1]);
        X(:,:,tr) = tmp;
    end

    xSZ = size(X);
    X = reshape(X,[xSZ(1) prod(xSZ(2:3))]);
    
    [S C U E L ERR LAM] = PCA_FIT_FULL_T(X,15);
    X = C;
    X = reshape(X,[15 xSZ(2:3)]);

    X = reshape(X,[588 xSZ(2:3)]);




    N = sum(X.*X,2);
    X = bsxfun(@times,X,N.^-1);

    %{
    options = optimset('Display','iter','PlotFcns',@optimplotfval);
    h = figure;
    x0 = ones(size(X,1),1)*size(X,1)^-1;
    func = @(K)myNewConvRegression2(K,X,Y',SZ,size(NEWX),P,h);
    x = fminsearch(func,x0,options);
    %}

    options = optimset('Display','iter','PlotFcns',@optimplotfval,'DiffMinChange',.1);
    x0 = ones(size(X,1),1)*size(X,1)^-1;
    h = figure;
    func = @(K)myNewConvRegression2(K,X,Y',SZ,size(NEWX),P,h,E,U);
    x = fmincon(func,x0,[],[],[],[],-ones(size(x0)),ones(size(x0)),@(q)myconR2(q),options);
%}

