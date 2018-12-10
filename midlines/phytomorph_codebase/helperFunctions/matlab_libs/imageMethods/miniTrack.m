function [T] = miniTrack(I,para)
        try 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % I(:,:,1)      is the first image
            % I(:,:,2)      is the second image
            % para.P        is the point to track
            % para.domain   is the tracking domain
            % para.RADIUS   is the clipping window
            % para.init_T   is the init vector for transformation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % clip out the image patch around the point
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [U1 U2] = ndgrid(para.P(1)-para.RADIUS:para.P(1)+para.RADIUS,para.P(2)-para.RADIUS:para.P(2)+para.RADIUS);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % external vs internal
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if para.externalInterp
                I(:,:,1) = ba_interp2(I(:,:,1),U2,U1);
                I(:,:,2) = ba_interp2(I(:,:,2),U2,U1);
            else
                I(:,:,1) = interp2(I(:,:,1),U2,U1);
                I(:,:,2) = interp2(I(:,:,2),U2,U1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % take the gradient
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [D1 D2] = gradient(I(:,:,1));
            dX = [1 1];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interpolate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if para.externalInterp
                Ii = ba_interp2(I(:,:,1),para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
                D1i = ba_interp2(D1,para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
                D2i = ba_interp2(D2,para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
            else
                Ii = interp2(I(:,:,1),para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
                D1i = interp2(D1,para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
                D2i = interp2(D2,para.domain(:,1)+para.RADIUS+1,para.domain(:,2)+para.RADIUS+1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % interpolate
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            D = [D1i D2i];

            icnt = 1;
            flag = 1;
            N = [];

            while flag & norm(dX) > para.threshold
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % init transformation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                TR = reshape(para.init_T,[2 3]);
                Xt = (TR*[para.domain ones(size(para.domain,1),1)]')';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % if internal
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if para.externalInterp
                    Gi = ba_interp2(I(:,:,2),Xt(:,1)+para.RADIUS+1,Xt(:,2)+para.RADIUS+1);
                else
                    Gi = interp2(I(:,:,2),Xt(:,1)+para.RADIUS+1,Xt(:,2)+para.RADIUS+1);    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % solution vector
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Mi = [D.*repmat(para.domain(:,1),[1 2]) D(:,1) D.*repmat(para.domain(:,2),[1 2]) D(:,2)];
                dY = Mi\(Ii-Gi);
                dX = [dY(1) dY(2) dY(4) dY(5) dY(3) dY(6)]';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % displace vector
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                para.init_T = para.init_T + dX';
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % measure image distance
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                N(icnt) = norm(Ii(:)-Gi(:));
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % ensure that the norm is minimizing.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                if icnt >= 2
                    flag = N(icnt) <  N(icnt-1);
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                icnt = icnt + 1;
            end
        catch ME
            
        end
end
