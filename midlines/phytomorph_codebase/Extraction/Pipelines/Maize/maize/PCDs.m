function [distance] = PCDs(X,U,E,targetCurve,percent,thresh,spline)
    sz = size(targetCurve);
    for e = 1:size(X,1)
        %sourceCurve = PCA_BKPROJ(X(e,:),E,U);
        spline.coefs = reshape(X(e,:),size(spline.coefs));
        sourceCurve = fnval(spline,linspace(0,1,size(targetCurve,2)));
        %sourceCurve = reshape(sourceCurve,sz);
        [d tE] = distanceMetric(sourceCurve,targetCurve);
        fidx = find(tE > .85);
        
        %distance(e) = mean(d();
        
        distance(e) = sum(d < thresh & tE > .85)/numel(d);
        distance(e) = 1 - distance(e);
        %{
        rmidx = d > thresh;
        d(rmidx) = [];
        distance = size(targetCurve,2) - numel(d);
        %}
        %{
        d = sort(d);
        d = d(1:round(percent*numel(d)));
        rmidx = d > thresh;
        d(rmidx) = [];
        if numel(d) < 30
             distance(e) = 100;
        else
             distance(e) = max(d);
        end
        %}
    end
    [j,midx] = min(distance);
    %sourceCurve = PCA_BKPROJ(X(midx,:),E,U);
    %sourceCurve = reshape(sourceCurve,sz);
    spline.coefs = reshape(X(midx,:),size(spline.coefs));
    sourceCurve = fnval(spline,linspace(0,1,size(targetCurve,2)));
    plot(targetCurve(1,:),targetCurve(2,:),'r');
    hold on;
    plot(sourceCurve(1,:),sourceCurve(2,:),'b');
    drawnow
end

function [d tE] = distanceMetric(source,target)
    sN = closedCurve.getNormals(source,5);
    tN = closedCurve.getNormals(target,5);
    D = sum(sN.*tN,1);
    for e = 1:size(target,2)
        dist = bsxfun(@minus,target(:,e),source);
        dist = sum(dist.*dist,1).^.5;
        [d(e) sidx] = min(dist);
        tE(e) = D(sidx);
    end
end