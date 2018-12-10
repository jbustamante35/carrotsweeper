function [out] = morphometrics_analysis(M,para)
    % para.K_window := cwt curvature 
    % para.SNIP     := snip for common domain
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    out.K = [];     % kurvature
    out.T = [];     % tip angle
    out.B = [];     % basis frames
    out.L = [];     % length
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % loop over frames
    for tm = 1:size(M,3)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % create curve
        iM = arcLength(M(:,:,tm),'arcLen',0);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % construct domain
        if tm == 1
            len = size(iM,1) - para.SNIP;
        end
        % snip out the common domain for the tm-th midline
        iMc = iM(para.SNIP:len,:);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % loop over points on common domain
        radius = para.radius;
        tmpB = zeros(size(iMc,1),2,2);
        for sp = 1:size(iMc,1)
            tmpB(sp,:,:) = igetFrame_atP(iM,sp,radius);
        end
        tmpB =  shiftdim(tmpB,-1);
        out.B = cat(1,out.B,tmpB);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calc curvature
        para_k{1} = para.K_window;
        % curvature
        o = cwtK(iMc,para_k);
        % stack curvature
        out.K = cat(1,out.K,o.K);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % tip angle
        
        %pidx = para.pidx;
        
        %out.T(tm) = atan2(b(2,1),b(1,1));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % length
        out.L(tm) = size(iM,1);
        fprintf(['done@' num2str(tm) ':' num2str(size(M,3)) '\n']);
    end
end