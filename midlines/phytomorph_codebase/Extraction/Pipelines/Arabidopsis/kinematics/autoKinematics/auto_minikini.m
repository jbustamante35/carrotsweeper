function [pointList] = auto_minikini(FileList,para)
    try 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 
        pointList = para.pointList;
        domainPara = para.domainPara;
        TH = para.THRESH;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % frame vector
        TIME = para.TIME;
        if isempty(para.TIME);TIME = 1:(numel(FileList));end
        
        %{
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % set figure handle
        if isfield(para,'h')
            if ~isempty(para.h)
                h = para.h;
            else
                h = figure;
            end
        else
            h = figure;
        end
        %}
        h = para.h;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % init tracking vars - only allow for RAD via user            
        RAD = para.domainPara{1}.value{1}(2);
        BUFFER = 60;
        RADIUS = RAD + BUFFER;
        THRESH = 10^-6;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create track domain
        domainPara = genDomains(para.domainPara);
        CONST = domainPara{1}.d(1:2,:)';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % for each frame
        for t = 1:(numel(TIME)-1)
            tm = TIME(t);
            tm_next = TIME(t+1);
            % init the track for non-linear
            initP = [1 0 0 1 0 0];
            % get the images from the tensor via the default view
            I = myReader(FileList{tm},'toGray',1);
            G = myReader(FileList{tm_next},'toGray',1);
            toTrack = pointList(:,:,t);
            tracked = zeros(size(toTrack));
            % track each point
            parfor pt = 1:size(pointList,1)
                tm = clock;
                % conditional statement for tracking
                if toTrack(pt,2) > TH
                    T = TR(I,G,toTrack(pt,:),RADIUS,CONST,THRESH,initP,1); % track is not tensor ready
                else
                    T = [0 0 0 0 0 0];
                end
                etm = etime(clock,tm);
                tracked(pt,:) = toTrack(pt,:) + fliplr((T(end-1:end)));
                fprintf(['tracked:' num2str(pt) ':' num2str(size(pointList,1)) ':' num2str(t) ':' num2str(numel(TIME)-1) ':' num2str(etm) '\n']);
            end
            pointList(:,:,t+1) = tracked;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % display to screen
            if ~isempty(h)                
                figure(h);
                imshow(I,[])            
                hold on
                scatter(pointList(:,2,t),pointList(:,1,t),'r.');
                for k = 1:size(pointList,1)
                    if pointList(k,2,t) > TH
                        plot(pointList(k,2,t) + RAD*cos(linspace(-pi,pi,200)),pointList(k,1,t) + RAD*sin(linspace(-pi,pi,200)),'r');
                    else
                        plot(pointList(k,2,t) + RAD*cos(linspace(-pi,pi,200)),pointList(k,1,t) + RAD*sin(linspace(-pi,pi,200)),'b');
                    end
                end
                hold off
                drawnow
            end
            % display to screen
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ~isempty(h);close(h);end
    catch ME
        ME.message
        ME.stack
        if ~isempty(h);close(h);end
    end
end

%%%%%%%%%%%%%%%%%%%%
% perform the tracking
function [T] = TR(I,G,P,RADIUS,X,PER,T,externalInterp)
    %%%%%%%%%%%%%%%%%%%%
    % this function is NOT tensor ready
    % it will drop back to OLD way
    % I is the first image
    % G is the second image
    % P is the point
    % RADIUS is the clipping window    
    %%%%%%%%%%%%%%%%%%%%
    % clip out the image patch around the point
    [U1 U2] = ndgrid(P(1)-RADIUS:P(1)+RADIUS,P(2)-RADIUS:P(2)+RADIUS);
    %%%%%%%%%%%%%%%%%%%%
    % external vs internal
    if externalInterp
        I = ba_interp2(I,U2,U1);
        G = ba_interp2(G,U2,U1);
    else
        I = interp2(I,U2,U1);
        G = interp2(G,U2,U1);
    end
    %%%%%%%%%%%%%%%%%%%%
    % take the gradient
    [D1 D2] = gradient(I);
    dX = [1 1];
    %%%%%%%%%%%%%%%%%%%%
    % interpolate
    if externalInterp
        Ii = ba_interp2(I,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
        D1i = ba_interp2(D1,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
        D2i = ba_interp2(D2,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
    else
        Ii = interp2(I,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
        D1i = interp2(D1,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
        D2i = interp2(D2,X(:,1)+RADIUS+1,X(:,2)+RADIUS+1);
    end
    %%%%%%%%%%%%%%%%%%%%
    % interpolate
    D = [D1i D2i];
    icnt = 1;
    flag = 1;
    N = [];
    while flag & norm(dX) > PER
        
        TR = reshape(T,[2 3]);        
        Xt = (TR*[X ones(size(X,1),1)]')';        
        % if internal
        if externalInterp
            Gi = ba_interp2(G,Xt(:,1)+RADIUS+1,Xt(:,2)+RADIUS+1);
        else
            Gi = interp2(G,Xt(:,1)+RADIUS+1,Xt(:,2)+RADIUS+1);    
        end
        % solution vector
        Mi = [D.*repmat(X(:,1),[1 2]) D(:,1) D.*repmat(X(:,2),[1 2]) D(:,2)];
        dY = Mi\(Ii-Gi);
        dX = [dY(1) dY(2) dY(4) dY(5) dY(3) dY(6)]';        
        % displace vector
        T = T + dX';
        N(icnt) = norm(Ii(:)-Gi(:));        
        % ensure that the norm is minimizing.
        if icnt >= 2
            flag = N(icnt) <  N(icnt-1);
        end        
        icnt = icnt + 1;
    end    
end
