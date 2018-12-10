function [E,U] = generateCurveSegments(dB,S,sm)
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    About:      
                generateCurveSegments.m  (Inputs are relative to 1200dpi)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Dependency: 
                cwtK_closed_imfilter.m, circshift.m
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Variable Definition:
                dB:      The information is needed. 
                S:
                sm:
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %}
    try
        Nsegs = [];
        cnt = 1;
        for e = 1:numel(dB)
            % generate the normal and tangent space for the curve
            E = getNormalsAndTangent(dB{e}',sm);        
            % generate curve segments
            segs = genS(dB{e}',S,E);        
            Nsegs = [Nsegs squeeze(segs(:,2,:))];
            e;
        end  
        [SIM C U E L ERR LAM] = PCA_FIT_FULL(Nsegs',2);
    catch ME
        close all;
        getReport(ME)
        fprintf(['******error in:generateCurveSegments.m******\n']);
    end
end




% generate the normal and tangent space
function [segs] = genS(segment,segmentSize,E)
    halfBlock = (segmentSize-1)/2;    
    sz = size(segment);
    J = [segment';segment';segment'];
    tmp1 = im2col(J(:,1),[segmentSize 1]);
    tmp2 = im2col(J(:,2),[segmentSize 1]);    
    tmp1 = tmp1(:,sz(2)+1-halfBlock:sz(2)+sz(2)-halfBlock);
    tmp2 = tmp2(:,sz(2)+1-halfBlock:sz(2)+sz(2)-halfBlock);
    segs = cat(3,tmp1,tmp2);
    segs = permute(segs,[2 1 3]);
    segs = permute(segs,[2 3 1]);
    for e = 1:size(segs,3)
        sz = size(segs,1);
        U = segs((sz-1)/2,:,e);
        segs(:,:,e) = bsxfun(@minus,segs(:,:,e),U);
        segs(:,:,e) = (E(:,:,e)'*segs(:,:,e)')';
    end
end


% generate the normal and tangent space
function [E] = getNormalsAndTangent(segment,S)
    sz = size(segment);
    J = [segment';segment';segment'];
    % calculate curvature
    d1X1 = cwt(J(:,1),S,'gaus1');
    d1X2 = cwt(J(:,2),S,'gaus1');            
    T = cat(3,d1X1,d1X2);
    L = sum(T.*T,3).^.5;
    T = bsxfun(@times,T,L.^-1);
    N = cat(3,T(:,:,2),-T(:,:,1));
    N = squeeze(N)';
    N = N(:,sz(2)+1:sz(2)+sz(2));
    T = squeeze(T)';
    T = T(:,sz(2)+1:sz(2)+sz(2));
    E = cat(3,permute(T,[2 1]),permute(N,[2 1]));
    E = permute(E,[2 3 1]);
end
