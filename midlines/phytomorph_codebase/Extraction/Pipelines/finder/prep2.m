function [vec] = prep2(movie,blockSize,reSize)
    vec = zeros(size(movie,4),prod(blockSize)*size(movie,3));
    CH = size(movie,3);
    parfor f = 1:size(movie,4)
        image = movie(:,:,:,f);
        
        % for each re-size request
        for re = 1:numel(reSize)
            rImage = (imresize(image,reSize(re)));
            M = zeros(prod(blockSize),CH);
            % for each color channel
            for k = 1:size(rImage,3)
                rSTOP = size(rImage,1) - 2*(blockSize(1) - 1);
                cSTOP = size(rImage,2) - 2*(blockSize(2) - 1);
                %[n1 n2] = ndgrid(1:rSTOP,1:cSTOP);
                %GR = [n1(:) n2(:)];
                tImage = rImage(:,:,k);
                tImage = tImage - imfilter(tImage,fspecial('average',blockSize));
                S = zeros(blockSize);
                for c = 1:blockSize(2)
                    for r = 1:blockSize(1)
                        %tGR = bsxfun(@plus,GR,[r-1 c-1]);
                        %g = sub2ind(size(tImage),tGR(:,1),tGR(:,2));
                        %vec(f,cnt) = vec(f,cnt) + sum(tImage(g));
                        %tmp = rImage(r:(rSTOP+(r-1)),c:(cSTOP+(c-1)),k);
                        tmp = tImage(r:(rSTOP+(r-1)),c:(cSTOP+(c-1)));
                        %vec(f,cnt) = vec(f,cnt) + sum(sum(rImage(r:(rSTOP+(r-1)),c:(cSTOP+(c-1)),k),1),2);
                        S(r,c) = sum(tmp(:));
                    end
                end
                M(:,k) = S(:);
            end
            vec(f,:) = vec(f,:) + M(:)';
        end
        fprintf(['done with frame:' num2str(f) '\n']);
    end
end
