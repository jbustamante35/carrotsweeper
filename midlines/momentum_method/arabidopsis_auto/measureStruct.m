function [S] = measureStruct(S,NP,NPK,SNIP,sz)
    try
        [angle length K snip frameOut snipC] = measurePhenoTypesOverStack(S.frameData,NP,NPK,3,SNIP);
        S.angle = angle;
        S.length = length;
        S.K = -permute(K,[2 3 1]);
        S.SNIP = permute(snip,[2 4 1 3]);
        S.frame = permute(frameOut,[1 2 4 3]);
        S.snipC = permute(snipC,[2 4 1 3]);
        S.ME = [];
    catch ME
        fprintf(['error in measure struct']);
        ME.getReport()
        
        
        S.ME = ME;
        S.angle = zeros(sz);
        S.length = zeros(sz);
        S.K = [];
        S.SNIP = [];
        S.frame = [];
        S.snipC = [];
    end
end


%%%%%%%%%%%%%%%%%%%%%%%
function [angle,length,K,snipOut frameOut snipCout] = measurePhenoTypesOverStack(S,NP,NPK,WAVELET,SNIP)
    % S : struct
    % NP : number points for tip angle measurement
    % NPK : number of points to measure curvature over
    % WAVELET : wavelet number for der
    % SNIP : length of SNIP to extract
    for t = 1:numel(S)
        angle(t,:) = measureAngleOverObjects(S{t},NP);
        length(t,:) = measureLengthOverObjects(S{t});
        K(:,:,t) = measureKurvatureOverObjects(S{t},WAVELET,NPK);
        [snipOut(:,:,:,t) frameOut(:,:,:,t) snipCout(:,:,:,t)] = extractSnipOverObjects(S{t},SNIP);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%
function [snip E snipC] = extractSnipOverObjects(S,NP)   
    for e = 1:numel(S.midlines)
        [snip(:,:,e) E(:,:,e) snipC(:,:,e)] = extractSnipFromCurve(S.midlines(e).data,NP);
    end
end 
function [snip E C] = extractSnipFromCurve(curve,NP)
    snip = curve(:,1:NP);
    [E C] = extractBasisFromSnip(snip);
end
function [E C] = extractBasisFromSnip(snip)
    stable = 10;
    [S C U E L ERR LAM] = PCA_FIT_FULL(snip(:,1:stable)',2);
    
    U = snip(:,1)';
    dC = diff(snip(:,1:stable),1,2);
    dC = mean(dC,2);
    if sign(E(:,1)'*dC) < 1
        E(:,1) = -E(:,1);
    end
    E(1,2) = -E(2,1);
    E(2,2) = E(1,1);
    C = PCA_REPROJ(snip',E,U);
    E = [E U'];
    C = C';

    %{    
    [S C U E L ERR LAM] = PCA_FIT_FULL(snip',2);
    
    
    dC = diff(snip,1,2);
    dC = mean(dC,2);
    if sign(E(:,1)'*dC) < 1
        E(:,1) = -E(:,1);
    end
    E(1,2) = -E(2,1);
    E(2,2) = E(1,1);
    C = PCA_REPROJ(snip',E,U);
    E = [E U'];
    C = C';
    %}
end
%%%%%%%%%%%%%%%%%%%%%%%
function [length] = measureLengthOverObjects(S)
    for e = 1:numel(S.midlines)
        length(e) = measureLengthFromCurve(S.midlines(e).data);
    end
end
function [length] = measureLengthFromCurve(curve)
    dC = diff(curve,1,2);
    dC = sum(dC.*dC,1).^.5;
    length = sum([0 dC]);
end
%%%%%%%%%%%%%%%%%%%%%%%
function [angle] = measureAngleOverObjects(S,NP)
    for e = 1:numel(S.midlines)
        angle(e) = measureAngleFromCurve(S.midlines(e).data,NP);
    end
end
function [angle] = measureAngleFromCurve(curve,NP)
    dC = diff(curve,1,2);
    dC = -mean(dC(:,1:NP),2);
    angle = atan2(dC(2),dC(1));
end
%%%%%%%%%%%%%%%%%%%%%%%
function [K] = measureKurvatureOverObjects(S,SMOOTH,NP)
    for e = 1:numel(S.midlines)
        K(e,:) = measureKurvatureFromCurve(S.midlines(e).data,SMOOTH,NP);
    end
end
function [K] = measureKurvatureFromCurve(curve,SMOOTH,NP)
    o = cwtK(curve',{SMOOTH});
    K = o.K(1:NP);
end