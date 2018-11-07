function [curve model labels] = extractCurves(img,model,init,final)
    SNIP = 15;
    modelSIZE = 21;
    if nargin < 3;init = 1;end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate doamin
    para.value{2} = [-pi pi 70];
    para.value{1} = [0 21 21];
    para.type = 'disk';    
    domain = phytoAdomain(para);
    domain.generateDomain();
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % handle flip
    img = handleFLIP(img,[]);
    % refine contour
    B = contourRefine(img,1000);
    % boundary traces
    dB = getContour(B);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create root label chain
    H = rootLabelProcessChain();
    contour = dB(1).contour;
    P = H(contour); 
    tipIndex = P{1}{2};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % create curve along boundary
    curve = phytoAcurve(fliplr(contour));
    % set default final
    if nargin == 4;if isempty(final);final = size(contour,1);end;else;final = size(contour,1);end
    % create curve along boundary
    curveSample = phytoAcurve(fliplr(contour(init:final,:)));
    % generate tangent space
    TM = curveSample.generateTangentSpace(10);
    % sample with domain
    sample = TM.sample(img,domain);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get data around tip for model
    if nargout == 2
        model = sample(tipIndex-modelSIZE:tipIndex+modelSIZE);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if model provided then isolate tip
    if nargin >= 2
        sample = sample(:);
        VEC = reshape(sample,[size(sample,1) size(sample,2)*size(sample,3)]);
        model = reshape(model,[size(model,1) size(model,2)*size(model,3)]);
        for e = 1:size(VEC,1)-size(model,1)
            tdM = VEC(e:e+size(model,1)-1,:) - model;
            dM(e) = sum(sum(tdM.*tdM,2),1);
        end
        [~,tidx] = min(dM);
    end
    if nargout == 3
        labels{1}{2} = tidx + init - 1 + (size(model,1)-1)/2;
        labels{1}{1} = curve(labels{1}{2});    
    end
    
    %{
    plot(contour(:,1),contour(:,2));
    hold on;
    plot(contour(tipIndex,1),contour(tipIndex,2),'r*')
    axis equal
    %}
    
end