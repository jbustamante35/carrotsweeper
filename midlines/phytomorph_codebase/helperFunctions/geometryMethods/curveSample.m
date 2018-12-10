function [G] = curveSample(imageName,curve,bundle,sampleParameter)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % this will sample a curve along the curve however, with
    % in the rotated frame given by bundle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUTS:   
    %           imageName           := file name for the image
    %           curve               := curve for sampling along
    %           bundle              := tangent and normal vector stack [t,n];
    %           sampleParameters    := sample parameters for the disk along the
    %                                   curve
    %                                    sampleParameters{1}.type = 'disk';
    %                                    sampleParameters{1}.value{1} = [rho_min rho_max num_rho];
    %                                    sampleParameters{1}.value{2} = [-pi pi num_phi];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUTS:  
    %           BV      : = a sequence of basis vectors for the tangent and
    %           normal space
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    disp = 1;
    
    % auto call to arc length parameterize curve
    curve = arcLength(curve,'arcLen',1);
    % number of points along curve
    numP = size(curve,1);
    % generate domains
    domains = genDomains(sampleParameter);
    % bump for extended read
    BUMP = 20;
    % return
    G = zeros([domains{1}.sz numP]);
    % loop over each point and sample
    for pt = 1:numP
        % create call to read image patch from image
        para{1} = curve(pt,:);
        para{2} = sampleParameter{1}.value{1}(2) + BUMP;
        para{3} = sampleParameter{1}.value{1}(2) + BUMP;        
        Is = myReader(imageName,para);
        
        tmp_frame = squeeze(bundle(pt,:,:));
        affine_transform = [tmp_frame' [para{2};para{3}]];
        tDomain = affine_transform*domains{1}.d;
        
        Y = myInterp(double(Is),tDomain');
        G(:,:,pt) = reshape(Y,domains{1}.sz);
        % if disp
        if disp
            imshow(G(:,:,pt),[])
            drawnow
        end 
    end
    if disp;close all;end
    G = permute(G,[3 1 2]);
end