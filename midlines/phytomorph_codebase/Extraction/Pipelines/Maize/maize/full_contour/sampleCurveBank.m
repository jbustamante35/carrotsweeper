function [atomicCurveSegment] = sampleCurveBank(I,curve,NH,sz,segmentSize,containerClass)    
    % sample along the curve
    [sam TAN] = sampleCurve(I,curve,NH);
    sam = reshape(sam,[sz(1) sz(2) size(curve.segment,2)]);
    
    TAN = circshift(TAN,[0 0 -(segmentSize-1)/2]);
    sam = circshift(sam,[0 0 -(segmentSize-1)/2]);
    
    tmp = [curve.segment];
    for e1 = 1:size(tmp,2)
        % get curve segment
        sig = tmp(:,1:segmentSize)';
        % fit with least squares spline
        fn = spap2(1,3,[1:size(sig,1)]',sig');
        % get function derivative
        fn2 = fnder(fn,1);
        % eval the der for tangengent vec
        tvec = fnval(fn2,(segmentSize-1)/2);
        tvec = tvec/norm(tvec);
        nvec = [tvec(2);-tvec(1)];
        % get the reference frame
        E = [tvec nvec];
        % get the center point
        U = fnval(fn,(segmentSize-1)/2);
        C = PCA_REPROJ(sig,E,U');            

        tmp = circshift(tmp,[0 -1]);

        atomicCurveSegment(e1) = containerClass();
        atomicCurveSegment(e1).orientation = TAN(:,:,e1);        
        atomicCurveSegment(e1).segment = sig';
        atomicCurveSegment(e1).Osegment = C';
        atomicCurveSegment(e1).image = sam(:,:,e1);
        atomicCurveSegment(e1).centerPoint = U;
    end
    atomicCurveSegment = circshift(atomicCurveSegment,[0 (segmentSize-1)/2]);
   
end

function [sam TAN] = sampleCurve(I,curve,NH)
    sam = zeros([size(NH,2) size(curve.segment,2)]);
    for i = 1:size(curve.segment,2)
        [sam(:,i) TAN(:,:,i)] = samplePoint(I,curve,i,NH);
    end
end

function [sam TAN] = samplePoint(I,curve,p,NH)
    E = getTangentSpace(curve,p);
    P = curve.segment(:,p);
    NH = E*NH;
    NH = bsxfun(@plus,NH,P);
    sam = ba_interp2(I,NH(1,:),NH(2,:));
    TAN = E;
end

% get tangent space
function [E] = getTangentSpace(curve,p)
    SNIP = 50;
    sig = curve.segment;
    sig = circshift(sig,[0 -(p-1)]);
    sig = sig(:,1:SNIP)';
    
    % curve fit method
    fn = spap2(1,3,[1:size(sig,1)]',sig');
    fn1 = fnder(fn,1);
    tvec = fnval(fn1,(SNIP-1)/2);
    tvec = tvec/norm(tvec);
    nvec = [tvec(2);-tvec(1)];
    E = [tvec nvec];
    
    %{
    % gradient method    
    E = gradient(curve.segment);
    E = E(:,p);
    E = [E [-E(2);E(1)]];
    %}
end

function [] = viewsam(I,curve)
    figure;
    for c = 1:numel(curve)
        imshow(I,[])
        for i = 1:size(curve(c).sam,3)
            plot(curve(c).segment(1,:),curve(c).segment(2,:))        
            hold on
            plot(curve(c).segment(1,i),curve(c).segment(2,i),'*')
            imshow(curve(c).sam(:,:,i),[])
            hold off
            drawnow
        end
    end
end

