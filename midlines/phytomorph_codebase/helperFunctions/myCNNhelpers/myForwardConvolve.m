function [Z] = myForwardConvolve(X, W,padTop, padLeft,padBottom, padRight,strideHeight, strideWidth,toFFT)
    if nargin == 8
        toFFT = 3;
    end
    % convolveForward2D   Convolve input images with kernels
    %
    % Inputs:
    % X - The input feature maps for a set of images. A (H)x(W)x(C)x(N) array.
    % W - The kernels for convolution. A (R)x(S)x(C)x(K) array.
    % padTop - Padding on the top.
    % padLeft - Padding on the left.
    % padBottom - Padding on the bottom.
    % padRight - Padding on the right.
    % strideHeight - The stride in the y direction.
    % strideWidth - The stride in the x direction.
    %
    % Output:
    % Z - The output feature maps for the images. A
    % floor((H + padTop + padBottom - R)/strideHeight + 1) x
    % floor((W + padLeft + padRight - S)/strideWidth + 1) x
    % (C) x (N) array.
    %
    % This corresponds to the cuDNN function "cudnnConvolutionForward".

    %   Copyright 2015-2017 The MathWorks, Inc.

    
    %toFFT = 0;
    
    imageHeight = size(X,1);
    imageWidth =  size(X,2);
    numInputMaps = size(X,3);
    numExamples =  size(X,4);

    % Apply the padding to the images if necessary
    if (padTop > 0) || (padLeft > 0) || (padBottom > 0) || (padRight > 0)
        if size(W,1) > 1
          %  X = iPadArray(X, padTop, 0, padBottom, 0);
        end
    end

    filterHeight = size(W,1);
    filterWidth = size(W,2);
    assert(size(W,3) == numInputMaps, 'Kernel dim 3 does not match data');
    numOutputMaps = size(W,4);

    %convolvedImageHeightWithoutStride = imageHeight + padTop + padBottom - filterHeight + 1;
    convolvedImageHeightWithoutStride = imageHeight;
    %convolvedImageWidthWithoutStride = imageWidth + padLeft + padRight - filterWidth + 1;
    convolvedImageWidthWithoutStride = imageWidth;

    % Allocate memory for the output
    Z = zeros(convolvedImageHeightWithoutStride, ...
        convolvedImageWidthWithoutStride, ...
        numOutputMaps, numExamples, 'like', X);

    % Flip the kernel in every dimension
    W = double(rot90(W,2));
    
    if toFFT == 1
        xsize = [size(X,1) size(X,2)];
        msize = [size(W,1) size(W,2)];
        for k = 1:numOutputMaps
            for c = 1:numInputMaps
                m = W(:,:,c,k);
                % pad m with zeros
                if any(msize < xsize)  % test, as user may have optimised already
                    m = exindex(m, 1:xsize(1), 1:xsize(2), {0});
                end
                % recentre m so that y(1,1) corresponds to mask centred on x(1,1)
                mc = 1 + floor(msize/2);
                me = mc + xsize - 1;
                mW(:,:,c,k) = fft2(exindex(m, mc(1):me(1), mc(2):me(2), 'circular'));
            end
        end
    end
    
    
    
    %{
    % added line
    W = double(W);
    [R T] = ndgrid(linspace(0,size(X,1),size(X,1)),linspace(-pi,pi,size(X,2)));
    X1 = R.*cos(T)+size(X,1)/2;
    X2 = R.*cos(T)+size(X,1)/2;
    % added line
    %}
    %numExamples
    %tic
  %{
    sz = floor((size(W)-1)/2);
    if any(mod(sz(1:2),2)==0)
        
    end
    %}
    
    if toFFT == 1
        for n = 1:numExamples
            for c = 1:numInputMaps
                X(:,:,c,n) = fft2(X(:,:,c,n));
            end
        end
    end
    
    for n = 1:numExamples
        for k = 1:numOutputMaps
            for c = 1:numInputMaps
                
                
                
                %tmp = fft2(X(:,:,c,n));
                %y = conv_fft2(tmp, mW(:,:,c,k), 'wrap',1);
                %y2 = imfilter(X(:,:,c,n),W(:,:,c,k),'circular','conv');
                
                %Z(:,:,k,n) = Z(:,:,k,n) + conv_fft2(X(:,:,c,n),mW(:,:,c,k),'circular','conv');
                if toFFT == 1
                    Z(:,:,k,n) = Z(:,:,k,n) + conv_fft2(X(:,:,c,n),mW(:,:,c,k),'wrap',1);
                elseif toFFT == 0
                    Z(:,:,k,n) = Z(:,:,k,n) + imfilter(X(:,:,c,n),W(:,:,c,k),'circular','conv');
                elseif toFFT == 2
                    Z(:,:,k,n) = Z(:,:,k,n) + conv_fft2(X(:,:,c,n),W(:,:,c,k),'wrap',0);
                elseif toFFT == 3
                    Z(:,:,k,n) = Z(:,:,k,n) + conv2(abs(fft(X(:,:,c,n),[],2)),W(:,:,c,k),'same');
                end
                
                %
                %Z(:,:,k,n) = Z(:,:,k,n) + conv2(padarray(X(:,:,c,n),sz(1:2),'circular','both'),W(:,:,c,k),'valid');
                %{
                F = ba_interp2(double(X(:,:,c,n)),X1,X2);
                Z(:,:,k,n) = Z(:,:,k,n) + imfilter(F,W(:,:,c,k),'circular','conv');
                % removed line
                %}
                % Perform 3D convolution
                %Z(:,:,k,n) = Z(:,:,k,n) + conv2(X(:,:,c,n), W(:,:,c,k), 'valid');


            end
        end
    end
    %toc
    % Downsample the result to account for stride
    Z = Z(1:strideHeight:end,1:strideWidth:end, :, :);
end


function Y = iPadArray(X, padTop, padLeft, padBottom, padRight)
paddedSize = size(X);
paddedSize(1) = paddedSize(1) + padTop + padBottom;
paddedSize(2) = paddedSize(2) + padLeft + padRight;
Y = zeros(paddedSize, 'like', X);
imageTop = padTop + 1;
imageBottom = padTop + size(X,1);
imageLeft = padLeft + 1;
imageRight = padLeft + size(X,2);
Y(imageTop:imageBottom, imageLeft:imageRight, :, :) = X;
end