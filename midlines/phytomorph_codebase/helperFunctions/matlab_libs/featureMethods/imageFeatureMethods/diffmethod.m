function [D] = diffmethod(I,para)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % take gradient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % INPUT: 
    %           I               := image
    %           para.method     := method for difference    
    %           para.sz         := size for cwt if hilbert
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OUTPUT: 
    %           D =             := gradient stacked with image
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % NOTES: add any diff method that i want to use to this file to try out different ones    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % switch on gradient type
    switch para.method.value
        % via finite differences
        case 'finite'
            [D1 D2] = gradient(I);
        % via wavelet
        case 'hilbert'
            D1 = zeros(size(I));
            D2 = zeros(size(I));
            for i = 1:size(I,1)
                D1(i,:) = -cwt(I(i,:),para.sz.value,'gaus1');
            end
            for i = 1:size(I,2)
                D2(:,i) = -cwt(I(:,i),para.sz.value,'gaus1')';
            end
        % gauss ~= wavelet
        case 'gauss'
            h = fspecial('gaussian',round(2*para.sz.value), para.sz.value);
            I = imfilter(I,h,'replicate');
            [D1 D2] = gradient(I);
    end
    % stack results
    D = cat(3,D1,D2);
end