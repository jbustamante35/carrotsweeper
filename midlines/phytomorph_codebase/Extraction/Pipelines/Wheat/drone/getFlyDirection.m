function [] = getFlyDirection(FileList)

    


    I = imread(FileList{1});
    
    
   
    [I,BOX] = imcrop(I);
    
    dX(1) = size(I,2)/2;
    dX(2) = size(I,1)/2;
    RAD = 400;
    [n2 n1] = ndgrid(linspace(0,RAD,100),linspace(0,pi,100));
    X = n2.*cos(n1) + size(I,2)/2;
    Y = n2.*sin(n1) + size(I,1)/2;
    sheet = cat(3,X,Y);
    
    
    stack = grabStack(FileList,1,1,1+75,BOX);
  
    [F] = lineSample(stack,sheet);
    
    
    %stack = double(stack);
    
    for e = 27:400
        
        
        [stack] = pushPop(FileList,e,stack,1,BOX);
        F = lineSamplePop(stack,F,sheet);
        
       
        
        
        g = [];
        KS = [];
        for k = 1:size(F,2)
            tmp = F(:,k,:,:);
            tmp = permute(tmp,[1 4 3 2]);
            tmp = rgb2gray(tmp/255)';
            
            
            E = edge(tmp);
            [g1 g2] = gradient(tmp);
            A = atan2(-g2(find(E)),g1(find(E)));
            A = atan2(-g2(:),g1(:));
            A(A<-pi/2) = A(A<-pi/2) + pi;
            A(A>pi/2) = A(A>pi/2) - pi;
            [f,xi] = ksdensity(A,linspace(-pi,pi,1000));
            KS = [KS;f];
            [~,midx]= max(f);
            g(k) = mean(A);
            g(k) = xi(midx);
        end
        [~,m] = min(g);
        
        imshow(stack(:,:,:,round((end-1)/2))/255,[])
        hold on
        
        for s = 1:size(sheet,2)
            plot(sheet(:,s,1),sheet(:,s,2),'r')
        end
        plot(sheet(:,m,1),sheet(:,m,2),'c')
        plot(sheet(:,1,1),sheet(:,1,2),'g')
        Q2 = squeeze(sheet(:,m,1:2));
        Q2 = fliplr(bsxfun(@plus,Q2,-dX));
        Q2(:,1) = -Q2(:,1);
        plot(Q2(:,1)+dX(1),Q2(:,2)+dX(2),'b')
        
        
        
        %imshow(tmp/255,[]);
        drawnow
    end
end