function [sI] = sampleAlongBase(point,imageStack,Domain,frameIndex)
    
    % sample image at point P
    % if not spec then the third value of point
    % is the 
    if nargin == 3;frameIndex = 3;end
    
    % buffer around the point for cropping
    BUFFER = 10;
    
    % get the image name
    imgF = imageStack{point(frameIndex)};
    
    % convert data to a point
    baseP = phytoApoint([point(1:2) 1]);
    % convert point to affine
    affineT = phytoAaffine.contructFromApoint(baseP);
    % move the doamin
    nD = affineT*Domain;
    
    % make the crop box based on location
    m = round(min(nD.d,[],2) - BUFFER);
    M = round(max(nD.d,[],2) + BUFFER);
    CP = mean([m M],2);
    SZ = round((M - m)/2);
    
    % make odd crop box
    if mod(SZ(1),2) == 1
        SZ(1) = SZ(1) + 1;
    end
    % make odd crop box
    if mod(SZ(2),2) == 1
        SZ(2) = SZ(2) + 1;
    end
    % read aaround crop box
    
    para{1} = round((CP(1:2)));
    para{2} = round(SZ(1));    
    para{3} = round(SZ(2));
    I = myReader(imgF,'atP',para);

    
    delta = (size(I)+1/2) + delta;
    
    
end
