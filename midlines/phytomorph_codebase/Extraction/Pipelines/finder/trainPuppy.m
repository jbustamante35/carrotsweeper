function [d] = trainPuppy(x,movie,target)
    tic;
    fprintf(['start with puppy\n']);
    x = reshape(x,[3 21*21]);

    reSize = linspace(.1,.4,5);
    gridSize = [21 21 3];
   
    
    func = @(m)myOp1(m,x);
    
    
    for fr = 1:size(movie,4)
        tic;
        d(fr) = conI(movie(:,:,:,fr),reSize,gridSize(1:2),func);
        toc
        %fprintf(['done with frame in:' num2str(toc) '\n']);
    end
    
    
    d = d - target';
    d = norm(d);
    fprintf(['done with puppy:' num2str(toc) '\n']);
end


