function [Tout Y] = pbPCA(I,T,X,meth)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % pbPCA - patch balanced PCA
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           I       : = image for interpolation
    %           T       : = transformation for the domain X
    %           X       : = domain
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:
    %           Tout    = reference frame made up by the eigenvectors
    %           Y       = sampling in the reference frame Tout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    para{1} = T;
    para{2} = X;
    R = [[cos(pi/2) sin(pi/2) 0];[-sin(pi/2) cos(pi/2) 0];[0 0 1]];
    %nX = (T*R*X)';
    nX = (T*X)';
    % added so that I can be file name
    Y = myReader(I,'iatP',para,'toGray',1);
    %I = myReader(I,'toGray');
    % interpolate
    %Y = myInterp(I,nX);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % threshold and pca on the disk
    Y = bindVec(Y);
    level = graythresh(Y(:));    
    idx = find(meth(Y,level));
    subX = nX(idx,:);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % call for PCA_FIT_FULL
    [S C U BV L ERR LAM] = PCA_FIT_FULL(subX,2);
    
   
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    BV(:,2) = twistVec(BV(:,1));
    Tout = [BV U'];
    delta = normalizeVec((T(:,end) - U'));
    targetFrame = [delta twistVec(delta) U'];
    Tout = minimizeTwist(Tout,targetFrame);
    
    %Tout(:,1:2) = inv(R)*Tout(:,1:2);
    %{
    figure;plot(subX(:,1),subX(:,2),'r.');hold on
    quiver(U(1),U(2),Tout(1),Tout(2),30)
    axis equal
    drawnow
    pause(1);
    close all;
    %}
    
    
    if nargout == 2
        para{1} = Tout;
        para{2} = X;    
        % added so that I can be file name
        Y = myReader(I,'iatP',para,'toGray');
    end
end

