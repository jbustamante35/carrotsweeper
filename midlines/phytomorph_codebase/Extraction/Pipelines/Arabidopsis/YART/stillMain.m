function [F] = stillMain(im1,im2,disp)
    % find corners
    [cim, r1, c1] = harris(im1,5, .05, 20);
    [cim, r2, c2] = harris(im2,5, .05, 20);
    % match via corr
    w = 21;    % Window size for correlation matching
    [m1,m2] = matchbycorrelation(im1, [r1';c1'], im2, [r2';c2'], w);
    % make homog    
    x1 = [m1(2,:); m1(1,:); ones(1,length(m1))];
    x2 = [m2(2,:); m2(1,:); ones(1,length(m1))]; 
    
    
    if disp
        close all
        % Display both images overlayed with inlying matched feature points
        imshow(double(im1)+double(im2),[]);hold on
        plot(m1(2,:),m1(1,:),'ro');
        plot(m2(2,:),m2(1,:),'go');
        for n = 1:size(x1,2)
            line([x1(1,n) x2(1,n)], [x1(2,n) x2(2,n)],'color',[1 0 0])
        end
        drawnow    
        %waitforbuttonpress
    end
    
    
    
    
    t = .1;  % Distance threshold for deciding outliers    
    [F, inliers] = myFit(x1, x2, t);
    
    if disp
        %%
        close all
        % Display both images overlayed with inlying matched feature points
        imshow(cat(3,double(im1)/255,double(im2)/255,zeros(size(im1))),[]);hold on
  
        plot(m1(2,inliers),m1(1,inliers),'go');
        plot(m2(2,inliers),m2(1,inliers),'go');    
        
        for n = inliers
            line([m1(2,n) m2(2,n)], [m1(1,n) m2(1,n)],'color',[0 1 0])
        end


        pt = x1(:,inliers);
        ptN = transformPointsInverse(F, pt(1:2,:)')';
        %ptN = F*pt;


        %imshow(double(im1)+double(im2),[]);hold on
        plot(pt(1,:),pt(2,:),'r*')
        plot(ptN(1,:),ptN(2,:),'r*')        
        for n = 1:size(pt,2)
            line([pt(1,n) ptN(1,n)], [pt(2,n) ptN(2,n)],'color',[1 0 0])
        end
        drawnow
        %waitforbuttonpress
        
    end
    %F(3,:) = [0 0 1];
    %F = maketform('affine',F');
end