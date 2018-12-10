function [m] = grapeCluster_MaskMeasure(mask)
    m = table;
    mask = logical(mask);
    mask = bwareaopen(mask,200);
    mask = imfill(mask,'holes');
    R = regionprops(mask,'EulerNumber','Area','Eccentricity','MajorAxisLength','MinorAxisLength','ConvexArea','EquivDiameter','Perimeter');
    f = fields(R);
    for e = 1:numel(f)
        m{1,['m_' f(e)]} = mean([R.(f{e})]);
        v = std([R.(f{e})],1,2);
        if isempty(v)
            v = NaN;
        end
        m{1,['s_' f(e)]} = v;
    end
    m{1,'numObjects'} = numel(R);
end

