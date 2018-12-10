function mosaicScaledImges=prepare_mosiac_elements(mosiacImages, mosaicHeigth, mosaicWidth, neededNOfFiles)
if nargin<4
    neededNOfFiles=Inf;
end

%% Read all images, and for each do the following:
% (1) resize image keeping aspect ratio to fit at least one of the mosaic dimentions 
% (2) crop to [mosaicHeigth,mosaicWidth] around resized image center
h = waitbar(0,'Reading mosaic images.');
nFiles=min( length(mosiacImages), neededNOfFiles);
mosaicScaledImges=cell(1,nFiles);
for iMosaic=1:nFiles
    img=mosiacImages{iMosaic};
    [img_rows, img_cols, clr_dim]=size(img);
    ratioH=mosaicHeigth/img_rows;
    ratioW=mosaicWidth/img_cols;
    resizeScale=max(ratioH, ratioW);
    res_img=uint8( imresize(img, resizeScale) );
    [resImgRows, resImgCols, clr_dim]=size(res_img);
    left=max(1, floor( (resImgCols-mosaicWidth)/2 ));
    top=max(1, floor( (resImgRows-mosaicHeigth)/2 ));
    cropRect=[left, top, mosaicWidth-1, mosaicHeigth-1];
    res_img=imcrop(res_img, cropRect);
%     if img_rows>img_cols
%         diff=img_rows-img_cols;
%         img=img(1+ceil(diff/2):end-floor(diff/2),:,:);  
%     elseif img_cols>img_rows
%         diff=img_cols-img_rows;
%         img=img(:,1+ceil(diff/2):end-floor(diff/2),:);  
%     end
%     res_img=uint8(imresize(img,[mosaicHeigth, mosaicWidth]));
    if clr_dim==1
            mosaicScaledImges{iMosaic}=repmat(res_img,[1,1,3]);
    else
            mosaicScaledImges{iMosaic}=res_img;
    end
    waitbar(iMosaic / nFiles,h);
end
close(h);