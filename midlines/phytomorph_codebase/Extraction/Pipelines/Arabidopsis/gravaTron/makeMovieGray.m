function [G] = makeMovieGray(I)
    for e = 1:size(I,4)
        G(:,:,e) = rgb2gray(double(I(:,:,:,e))/255);
        fprintf(['Done with gray for:' num2str(e) ':' num2str(size(I,4)) '\n']);
    end
end