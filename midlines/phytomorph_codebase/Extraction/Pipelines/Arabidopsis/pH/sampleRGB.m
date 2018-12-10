function [SAM]= sampleRGB(FileList,rootCurve,tip,SNIP,SMOOTH)
    for e = 1:numel(FileList)
        I = imread(FileList{e});
         I = flipdim(I,1);
        I = imfilter(I,fspecial('gaussian',[SMOOTH SMOOTH]));
        for c = 1:2
            SAM{1}(:,e,c) = ba_interp2(double(I(:,:,c)),rootCurve{e}(tip(e)-SNIP:tip(e),2),rootCurve{e}(tip(e)-SNIP:tip(e),1))/255;
            SAM{2}(:,e,c) = ba_interp2(double(I(:,:,c)),rootCurve{e}(tip(e):tip(e)+SNIP,2),rootCurve{e}(tip(e):tip(e)+SNIP,1))/255;
        end
    end
    SAM{1}(:,:,3) = SAM{1}(:,:,2).*SAM{1}(:,:,1).^-1;
    SAM{2}(:,:,3) = SAM{2}(:,:,2).*SAM{2}(:,:,1).^-1;
end