function [D] = myNewConvRegression(K,X,Y,SZ,P)
    for tr = 1:size(X,4)
        tmp = X(:,:,:,tr);
        tmp = im2colF(tmp,SZ,[1 1 1]);
        CONV = mtimesx(K',tmp);
        CONV = col2im(CONV,SZ(1:2),[size(X,1),size(X,2)]);
        Q(tr,1) = sum(sum(P(:,:,1)*CONV,1),2);
        Q(tr,2) = sum(sum(P(:,:,2)*CONV,1),2);
        
        
        if tr == 1
            imshow(CONV,[]);
            drawnow
        end
        
    end
    D = mean(sum((Q - Y).^2,2),1);
end

%{
    P = [];
    SZ = [7 7 3];
    HW = (SZ-1)/2;
    [P(:,:,1),P(:,:,2)] = ndgrid(HW(1)+1:(size(NEWX,1)-HW(1)),HW(2)+1:(size(NEWX,2)-HW(2)));
    
    options = optimset('Display','iter');

    x0 = ones(prod(SZ),1)*prod(SZ)^-1;



    for tr = 1:size(NEWX,4)
        tmp = NEWX(:,:,:,tr);
        tmp = im2colF(tmp,SZ,[1 1 1]);
        X(:,:,tr) = tmp;
    end

    func = @(K)myNewConvRegression(K,NEWX,Y',SZ,P);

    
    x = fminsearch(func,x0,options);

    x = fmincon(func,xO,[],[],ones(1,prod(SZ)),1);
%}

