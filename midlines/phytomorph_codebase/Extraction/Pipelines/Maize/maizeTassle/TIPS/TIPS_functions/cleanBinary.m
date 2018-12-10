function [ bin, thresh, perimeter] = cleanBinary( subtractedImage )
%CLEANBINARY binarizes the tassel image and removes all but the largest
%single object from the image (ideally, the tassel).
    gray = im2double(rgb2gray(subtractedImage));
    
    thresh = graythresh(gray);
    
    bin = gray > thresh;
    
    L = bwlabel(bin);
    R = regionprops(logical(bin));
    
    keep = find([R.Area] == max([R.Area]));
    bin(L ~= keep) = 0;
    
    bin = imclose(bin, strel('disk', 5)); % Close small holes first
    
    outline = bwmorph(bin, 'remove');
    L = bwlabel(outline);
    R = regionprops(outline);
    keep = find([R.Area] == max([R.Area]));
    outline(L ~= keep) = 0;

    
    perimeter = sum(sum(outline));

end

